program test_mqc_finite_differences
   !! Test finite difference geometry perturbation
   use pic_types, only: dp
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_finite_differences, only: displaced_geometry_t, generate_perturbed_geometries, &
                                     energy_hessian_task_count, energy_hessian_task, &
                                     energy_hessian_geometry, finite_diff_hessian_from_energies, &
                                     finite_diff_gradient_from_energies
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

   call check_energy_hessian(test_passed)

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
   end if

   print *, ""
   print *, "SOME TESTS FAILED"
   stop 1

contains

   subroutine check_energy_hessian(passed)
      !! The Hessian and gradient from energies, against a function with known derivatives
      !!
      !! `model_energy` is quadratic, plus three cubic terms and two quartic
      !! ones. The second differences are exact for the cubics, and for the
      !! quartics they err by `O(h**2)` in exactly two elements, `H(1,1)` and
      !! `H(1,2)`, which is checked by halving the step.
      logical, intent(inout) :: passed

      integer, parameter :: N_COORDS = 9
      real(dp), parameter :: STEP = 0.02_dp
      type(physical_fragment_t) :: reference, displaced
      real(dp), allocatable :: energies(:), hessian(:, :), half_hessian(:, :), gradient(:, :)
      real(dp) :: exact(N_COORDS, N_COORDS), exact_gradient(N_COORDS), x0(N_COORDS)
      real(dp) :: worst, error_full, error_half
      integer :: task, n_tasks, coord_i, sign_i, coord_j, sign_j, i, j
      integer :: seen_single, seen_pair, seen_reference

      call create_test_water_molecule(reference)
      reference%coordinates(:, 1) = [0.10_dp, -0.20_dp, 0.30_dp]
      reference%coordinates(:, 2) = [1.40_dp, 0.20_dp, -0.10_dp]
      reference%coordinates(:, 3) = [-0.30_dp, 1.60_dp, 0.25_dp]
      x0 = reshape(reference%coordinates, [N_COORDS])

      ! Task bookkeeping: 1 + 2 n**2 tasks, one reference, 2n singles, 4 per pair
      n_tasks = energy_hessian_task_count(N_COORDS)
      seen_reference = 0
      seen_single = 0
      seen_pair = 0
      do task = 1, n_tasks
         call energy_hessian_task(task, N_COORDS, coord_i, sign_i, coord_j, sign_j)
         if (coord_i == 0) then
            seen_reference = seen_reference + 1
         else if (coord_j == 0) then
            seen_single = seen_single + 1
         else if (coord_j > coord_i) then
            seen_pair = seen_pair + 1
         end if
      end do
      if (n_tasks /= 1 + 2*N_COORDS**2 .or. seen_reference /= 1 .or. &
          seen_single /= 2*N_COORDS .or. seen_pair /= 2*N_COORDS*(N_COORDS - 1)) then
         print *, "FAIL: energy task table has", n_tasks, "tasks:", seen_reference, seen_single, seen_pair
         passed = .false.
      else
         print *, "PASS: energy task table holds", n_tasks, "tasks for", N_COORDS, "coordinates"
      end if

      ! The geometry of every task is the one its decoded displacement names
      worst = 0.0_dp
      do task = 1, n_tasks
         call energy_hessian_task(task, N_COORDS, coord_i, sign_i, coord_j, sign_j)
         call energy_hessian_geometry(reference, task, STEP, displaced)
         worst = max(worst, maxval(abs(reshape(displaced%coordinates, [N_COORDS]) - x0 - &
                                       displacement_vector(coord_i, sign_i, coord_j, sign_j, STEP))))
         call displaced%destroy()
      end do
      if (worst > 1.0e-14_dp) then
         print *, "FAIL: a task geometry differs from its decoded displacement by", worst
         passed = .false.
      else
         print *, "PASS: every task geometry matches its decoded displacement"
      end if

      call exact_derivatives(x0, exact_gradient, exact)

      energies = energies_at(reference, STEP)
      call finite_diff_hessian_from_energies(N_COORDS, energies, STEP, hessian)
      energies = energies_at(reference, 0.5_dp*STEP)
      call finite_diff_hessian_from_energies(N_COORDS, energies, 0.5_dp*STEP, half_hessian)

      ! Symmetric, to the last bit
      worst = maxval(abs(hessian - transpose(hessian)))
      if (worst > 0.0_dp) then
         print *, "FAIL: Hessian from energies is not symmetric:", worst
         passed = .false.
      else
         print *, "PASS: Hessian from energies is symmetric"
      end if

      ! Away from H(1,1) and H(1,2) the quartic terms contribute nothing
      worst = 0.0_dp
      do j = 1, N_COORDS
         do i = 1, N_COORDS
            if (i == 1 .and. (j == 1 .or. j == 2)) cycle
            if (j == 1 .and. i == 2) cycle
            worst = max(worst, abs(hessian(i, j) - exact(i, j)))
         end do
      end do
      if (worst > 1.0e-8_dp) then
         print *, "FAIL: Hessian from energies off by", worst, "where the model is a cubic"
         passed = .false.
      else
         print *, "PASS: Hessian from energies exact to", worst, "for cubic terms"
      end if

      ! Where the quartics enter, the error shrinks as h**2
      error_full = max(abs(hessian(1, 1) - exact(1, 1)), abs(hessian(1, 2) - exact(1, 2)))
      error_half = max(abs(half_hessian(1, 1) - exact(1, 1)), abs(half_hessian(1, 2) - exact(1, 2)))
      if (error_full < 1.0e-6_dp .or. error_full > 1.0e-1_dp .or. &
          abs(error_full/error_half - 4.0_dp) > 0.2_dp) then
         print *, "FAIL: quartic error", error_full, error_half, "is not O(h**2)"
         passed = .false.
      else
         print *, "PASS: quartic error", error_full, "->", error_half, "falls as h**2"
      end if

      ! The gradient out of the same energies
      energies = energies_at(reference, STEP)
      call finite_diff_gradient_from_energies(reference%n_atoms, energies, STEP, gradient)
      worst = maxval(abs(reshape(gradient, [N_COORDS]) - exact_gradient))
      if (worst > 2.0e-3_dp) then
         print *, "FAIL: gradient from energies off by", worst
         passed = .false.
      else
         print *, "PASS: gradient from energies within", worst, "of the analytic one"
      end if

      call reference%destroy()
   end subroutine check_energy_hessian

   function energies_at(reference, step) result(energies)
      !! The model energy at every task geometry
      type(physical_fragment_t), intent(in) :: reference
      real(dp), intent(in) :: step
      real(dp), allocatable :: energies(:)

      type(physical_fragment_t) :: displaced
      integer :: task

      allocate (energies(energy_hessian_task_count(3*reference%n_atoms)))
      do task = 1, size(energies)
         call energy_hessian_geometry(reference, task, step, displaced)
         energies(task) = model_energy(reshape(displaced%coordinates, [3*reference%n_atoms]))
         call displaced%destroy()
      end do
   end function energies_at

   function displacement_vector(coord_i, sign_i, coord_j, sign_j, step) result(vector)
      !! The step a decoded task names, in the 9 flattened coordinates
      integer, intent(in) :: coord_i, sign_i, coord_j, sign_j
      real(dp), intent(in) :: step
      real(dp) :: vector(9)

      vector = 0.0_dp
      if (coord_i > 0) vector(coord_i) = sign_i*step
      if (coord_j > 0) vector(coord_j) = vector(coord_j) + sign_j*step
   end function displacement_vector

   pure function model_energy(x) result(energy)
      !! A quadratic form, three cubic terms and two quartic ones, in 9 variables
      real(dp), intent(in) :: x(9)
      real(dp) :: energy

      real(dp) :: a(9, 9), g(9)
      integer :: i, j

      call model_quadratic(a, g)
      energy = -76.0_dp + dot_product(g, x)
      do j = 1, 9
         do i = 1, 9
            energy = energy + 0.5_dp*a(i, j)*x(i)*x(j)
         end do
      end do
      energy = energy + 0.3_dp*x(1)**3 - 0.2_dp*x(1)*x(2)*x(4) + 0.5_dp*x(2)**2*x(5) &
               + 0.07_dp*x(1)**4 + 0.11_dp*x(1)**3*x(2)
   end function model_energy

   pure subroutine model_quadratic(a, g)
      !! The symmetric matrix and the linear term of `model_energy`
      real(dp), intent(out) :: a(9, 9), g(9)

      integer :: i, j

      do j = 1, 9
         do i = 1, 9
            a(i, j) = 0.05_dp/real(1 + abs(i - j), dp) + merge(0.4_dp + 0.03_dp*i, 0.0_dp, i == j)
         end do
         g(j) = 0.02_dp*real(j, dp) - 0.07_dp
      end do
   end subroutine model_quadratic

   subroutine exact_derivatives(x, gradient, hessian)
      !! The derivatives of `model_energy` at `x`, by hand
      real(dp), intent(in) :: x(9)
      real(dp), intent(out) :: gradient(9), hessian(9, 9)

      real(dp) :: a(9, 9), g(9)

      call model_quadratic(a, g)
      gradient = g + matmul(a, x)
      hessian = a

      gradient(1) = gradient(1) + 0.9_dp*x(1)**2 - 0.2_dp*x(2)*x(4) + 0.28_dp*x(1)**3 + 0.33_dp*x(1)**2*x(2)
      gradient(2) = gradient(2) - 0.2_dp*x(1)*x(4) + 1.0_dp*x(2)*x(5) + 0.11_dp*x(1)**3
      gradient(4) = gradient(4) - 0.2_dp*x(1)*x(2)
      gradient(5) = gradient(5) + 0.5_dp*x(2)**2

      hessian(1, 1) = hessian(1, 1) + 1.8_dp*x(1) + 0.84_dp*x(1)**2 + 0.66_dp*x(1)*x(2)
      hessian(1, 2) = hessian(1, 2) - 0.2_dp*x(4) + 0.33_dp*x(1)**2
      hessian(1, 4) = hessian(1, 4) - 0.2_dp*x(2)
      hessian(2, 4) = hessian(2, 4) - 0.2_dp*x(1)
      hessian(2, 2) = hessian(2, 2) + 1.0_dp*x(5)
      hessian(2, 5) = hessian(2, 5) + 1.0_dp*x(2)
      hessian(2, 1) = hessian(1, 2)
      hessian(4, 1) = hessian(1, 4)
      hessian(4, 2) = hessian(2, 4)
      hessian(5, 2) = hessian(2, 5)
   end subroutine exact_derivatives

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
