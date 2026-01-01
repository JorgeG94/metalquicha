program test_mqc_finite_differences
   !! Test finite difference geometry perturbation
   use pic_types, only: dp
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_finite_differences, only: displaced_geometry_t, generate_perturbed_geometries
   implicit none

   type(physical_fragment_t) :: ref_geom
   type(displaced_geometry_t), allocatable :: forward_geoms(:), backward_geoms(:)
   real(dp) :: displacement
   integer :: n_displacements, i
   logical :: test_passed

   print *, "Testing finite difference geometry generation..."

   ! Create a simple test geometry (H2O molecule)
   call create_test_water_molecule(ref_geom)

   ! Generate perturbed geometries
   displacement = 0.001_dp  ! Bohr
   call generate_perturbed_geometries(ref_geom, displacement, forward_geoms, backward_geoms)

   n_displacements = 3*ref_geom%n_atoms
   test_passed = .true.

   ! Verify we got the right number of geometries
   if (size(forward_geoms) /= n_displacements .or. size(backward_geoms) /= n_displacements) then
      print *, "FAIL: Expected", n_displacements, "displacements, got", &
         size(forward_geoms), "forward and", size(backward_geoms), "backward"
      test_passed = .false.
   else
      print *, "PASS: Generated", n_displacements, "forward and backward geometries"
   end if

   ! Verify displacement magnitudes and directions
   do i = 1, n_displacements
      if (abs(forward_geoms(i)%displacement - displacement) > 1.0e-12_dp) then
         print *, "FAIL: Forward displacement", i, "incorrect:", forward_geoms(i)%displacement
         test_passed = .false.
      end if
      if (abs(backward_geoms(i)%displacement - displacement) > 1.0e-12_dp) then
         print *, "FAIL: Backward displacement", i, "incorrect:", backward_geoms(i)%displacement
         test_passed = .false.
      end if
      if (forward_geoms(i)%direction /= 1) then
         print *, "FAIL: Forward direction", i, "incorrect:", forward_geoms(i)%direction
         test_passed = .false.
      end if
      if (backward_geoms(i)%direction /= -1) then
         print *, "FAIL: Backward direction", i, "incorrect:", backward_geoms(i)%direction
         test_passed = .false.
      end if
   end do

   if (test_passed) then
      print *, "PASS: All displacement magnitudes and directions correct"
   end if

   ! Verify that only the intended coordinate changed
   do i = 1, n_displacements
      block
         integer :: iatom, icoord, jatom, jcoord
         real(dp) :: expected_diff
         logical :: coord_test_passed

         iatom = forward_geoms(i)%atom_index
         icoord = forward_geoms(i)%coordinate

         coord_test_passed = .true.

         do jatom = 1, ref_geom%n_atoms
            do jcoord = 1, 3
               if (jatom == iatom .and. jcoord == icoord) then
                  ! This coordinate should be displaced
                  expected_diff = displacement
               else
                  ! All other coordinates should be unchanged
                  expected_diff = 0.0_dp
               end if

               if (abs(forward_geoms(i)%geometry%coordinates(jcoord, jatom) - &
                       ref_geom%coordinates(jcoord, jatom) - expected_diff) > 1.0e-12_dp) then
                  print *, "FAIL: Forward geom", i, "atom", jatom, "coord", jcoord, "incorrect"
                  coord_test_passed = .false.
               end if

               if (jatom == iatom .and. jcoord == icoord) then
                  expected_diff = -displacement
               else
                  expected_diff = 0.0_dp
               end if

               if (abs(backward_geoms(i)%geometry%coordinates(jcoord, jatom) - &
                       ref_geom%coordinates(jcoord, jatom) - expected_diff) > 1.0e-12_dp) then
                  print *, "FAIL: Backward geom", i, "atom", jatom, "coord", jcoord, "incorrect"
                  coord_test_passed = .false.
               end if
            end do
         end do

         if (.not. coord_test_passed) then
            test_passed = .false.
         end if
      end block
   end do

   if (test_passed) then
      print *, "PASS: All coordinate displacements are correct"
   end if

   ! Cleanup
   call ref_geom%destroy()
   do i = 1, n_displacements
      call forward_geoms(i)%destroy()
      call backward_geoms(i)%destroy()
   end do
   deallocate (forward_geoms, backward_geoms)

   if (test_passed) then
      print *, ""
      print *, "ALL TESTS PASSED"
      stop 0
   else
      print *, ""
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

      ! H atoms
      geom%element_numbers(2) = 1  ! Hydrogen
      geom%coordinates(:, 2) = [1.5_dp, 0.0_dp, 0.0_dp]

      geom%element_numbers(3) = 1  ! Hydrogen
      geom%coordinates(:, 3) = [0.0_dp, 1.5_dp, 0.0_dp]

   end subroutine create_test_water_molecule

end program test_mqc_finite_differences
