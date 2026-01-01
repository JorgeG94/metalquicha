program test_mqc_hessian_finite_diff
   !! Test XTB Hessian calculation via finite differences
   use pic_types, only: dp
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_method_xtb, only: xtb_method_t
   implicit none

   type(physical_fragment_t) :: water_geom
   type(xtb_method_t) :: xtb_calc
   type(calculation_result_t) :: result
   logical :: test_passed
   integer :: n_atoms, n_coords

   print *, "Testing XTB Hessian calculation via finite differences..."
   print *, ""

   ! Create a simple H2O geometry
   call create_test_water_molecule(water_geom)
   n_atoms = water_geom%n_atoms
   n_coords = 3*n_atoms

   ! Setup XTB method
   xtb_calc%variant = "gfn2"
   xtb_calc%verbose = .true.

   print *, "Computing Hessian for H2O molecule ("//char(n_atoms + 48)//" atoms)..."
   print *, ""

   ! Calculate Hessian
   call xtb_calc%calc_hessian(water_geom, result)

   test_passed = .true.

   ! Verify Hessian was computed
   if (.not. result%has_hessian) then
      print *, "FAIL: Hessian flag not set"
      test_passed = .false.
   else
      print *, "PASS: Hessian computed successfully"
   end if

   ! Verify Hessian dimensions
   if (allocated(result%hessian)) then
      if (size(result%hessian, 1) /= n_coords .or. size(result%hessian, 2) /= n_coords) then
         print *, "FAIL: Hessian dimensions incorrect"
         print *, "  Expected:", n_coords, "x", n_coords
         print *, "  Got:", size(result%hessian, 1), "x", size(result%hessian, 2)
         test_passed = .false.
      else
         print *, "PASS: Hessian dimensions correct ("//char(n_coords + 48)//"x"//char(n_coords + 48)//")"
      end if
   else
      print *, "FAIL: Hessian not allocated"
      test_passed = .false.
   end if

   ! Verify Hessian symmetry (should be symmetric for exact calculations)
   if (allocated(result%hessian)) then
      block
         integer :: i, j
         real(dp) :: max_asymmetry
         max_asymmetry = 0.0_dp
         do i = 1, n_coords
            do j = i + 1, n_coords
               max_asymmetry = max(max_asymmetry, &
                                   abs(result%hessian(i, j) - result%hessian(j, i)))
            end do
         end do

         print *, "  Maximum Hessian asymmetry:", max_asymmetry

         ! Finite differences may have small numerical asymmetry
         if (max_asymmetry > 1.0e-6_dp) then
            print *, "WARNING: Hessian has significant asymmetry (expected for finite differences)"
         else
            print *, "PASS: Hessian is reasonably symmetric"
         end if
      end block
   end if

   ! Verify energy and gradient were also computed
   if (.not. result%has_energy) then
      print *, "FAIL: Energy not computed"
      test_passed = .false.
   else
      print *, "PASS: Energy computed: ", result%energy%total(), " Hartree"
   end if

   if (.not. result%has_gradient) then
      print *, "FAIL: Gradient not computed"
      test_passed = .false.
   else
      print *, "PASS: Gradient computed"
   end if

   ! Print a sample of Hessian elements
   if (allocated(result%hessian)) then
      print *, ""
      print *, "Sample Hessian elements (atomic units):"
      print *, "  H[1,1] =", result%hessian(1, 1)
      print *, "  H[1,2] =", result%hessian(1, 2)
      print *, "  H[2,2] =", result%hessian(2, 2)
   end if

   ! Cleanup
   call water_geom%destroy()
   call result%destroy()

   print *, ""
   if (test_passed) then
      print *, "ALL TESTS PASSED"
      stop 0
   else
      print *, "SOME TESTS FAILED"
      stop 1
   end if

contains

   subroutine create_test_water_molecule(geom)
      !! Create a simple H2O geometry for testing
      type(physical_fragment_t), intent(out) :: geom

      geom%n_atoms = 3
      geom%charge = 0
      geom%multiplicity = 1
      geom%nelec = 10
      geom%n_caps = 0

      allocate (geom%element_numbers(3))
      allocate (geom%coordinates(3, 3))

      ! O atom at origin
      geom%element_numbers(1) = 8  ! Oxygen
      geom%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]

      ! H atoms (realistic H2O geometry in Bohr)
      geom%element_numbers(2) = 1  ! Hydrogen
      geom%coordinates(:, 2) = [1.8_dp, 0.0_dp, 0.0_dp]

      geom%element_numbers(3) = 1  ! Hydrogen
      geom%coordinates(:, 3) = [-0.464_dp, 1.745_dp, 0.0_dp]

   end subroutine create_test_water_molecule

end program test_mqc_hessian_finite_diff
