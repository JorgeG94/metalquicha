!! Hessians by finite differences of analytic gradients
module mqc_semi_numerical_hessian
   !! Builds a Hessian from central differences of a method's analytic
   !! gradients:
   !!
   !!   H[i,j] = (g_j(x_i + h) - g_j(x_i - h)) / 2h
   !!
   !! Semi-numerical: only one derivative is taken numerically, the other is
   !! analytic. Cost is 6N gradient evaluations for N atoms, plus the
   !! undisplaced point, each of them a full SCF.
   !!
   !! `finite_difference_hessian_energies` is the fully numerical form, second
   !! differences of energies alone, for a method with no analytic gradient.
   !! It costs `1 + 18 N**2` energy evaluations.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: ERROR_GENERIC, ERROR_VALIDATION
   use mqc_finite_differences, only: generate_perturbed_geometries, displaced_geometry_t, &
                                     finite_diff_hessian_from_gradients, &
                                     finite_diff_dipole_derivatives, DEFAULT_DISPLACEMENT, &
                                     energy_hessian_task_count, energy_hessian_task, &
                                     energy_hessian_geometry, finite_diff_hessian_from_energies, &
                                     finite_diff_gradient_from_energies
   implicit none
   private

   public :: finite_difference_hessian  !! Hessian from central differences of gradients
   public :: finite_difference_hessian_energies  !! Hessian from second differences of energies
   public :: suggest_numerical_hessian  !! Add the `numerical` route to a gradient refusal

contains

   subroutine finite_difference_hessian(method, fragment, result, verbose, displacement_in)
      !! Hessian of `method` at `fragment`, by central differences of gradients
      class(qc_method_t), intent(in) :: method            !! Any method with analytic gradients
      type(physical_fragment_t), intent(in) :: fragment   !! Reference geometry
      type(calculation_result_t), intent(out) :: result   !! Energy, gradient and Hessian
      logical, intent(in), optional :: verbose
      real(dp), intent(in), optional :: displacement_in
         !! Step in Bohr; `DEFAULT_DISPLACEMENT` when absent.

      type(displaced_geometry_t), allocatable :: forward_geoms(:), backward_geoms(:)
      real(dp), allocatable :: forward_gradients(:, :, :), backward_gradients(:, :, :)
      real(dp), allocatable :: forward_dipoles(:, :), backward_dipoles(:, :)
      type(calculation_result_t) :: point
      real(dp) :: displacement
      integer :: n_atoms, n_displacements, i
      logical :: is_verbose, have_dipoles

      is_verbose = .false.
      if (present(verbose)) is_verbose = verbose

      n_atoms = fragment%n_atoms
      n_displacements = 3*n_atoms
      displacement = DEFAULT_DISPLACEMENT
      if (present(displacement_in)) then
         if (displacement_in > 0.0_dp) displacement = displacement_in
      end if

      if (is_verbose) then
         call logger%info("Hessian by central differences of analytic gradients")
         call logger%info("  Atoms:                 "//to_char(n_atoms))
         call logger%info("  Gradient evaluations:  "//to_char(2*n_displacements))
         call logger%info("  Displacement:          "//to_char(displacement)//" Bohr")
      end if

      call generate_perturbed_geometries(fragment, displacement, forward_geoms, backward_geoms)
      allocate (forward_gradients(n_displacements, 3, n_atoms))
      allocate (backward_gradients(n_displacements, 3, n_atoms))
      ! The same displacements give dmu/dR, and so IR intensities.
      allocate (forward_dipoles(n_displacements, 3), backward_dipoles(n_displacements, 3))
      forward_dipoles = 0.0_dp
      backward_dipoles = 0.0_dp
      have_dipoles = .true.

      do i = 1, n_displacements
         call method%calc_gradient(forward_geoms(i)%geometry, point)
         if (point%has_error .or. .not. point%has_gradient) then
            call fail(result, point, "forward", i)
            return
         end if
         forward_gradients(i, :, :) = point%gradient
         if (point%has_dipole) then
            forward_dipoles(i, :) = point%dipole
         else
            have_dipoles = .false.
         end if
         call point%destroy()

         call method%calc_gradient(backward_geoms(i)%geometry, point)
         if (point%has_error .or. .not. point%has_gradient) then
            call fail(result, point, "backward", i)
            return
         end if
         backward_gradients(i, :, :) = point%gradient
         if (point%has_dipole) then
            backward_dipoles(i, :) = point%dipole
         else
            have_dipoles = .false.
         end if
         call point%destroy()
      end do

      call finite_diff_hessian_from_gradients(fragment, forward_gradients, backward_gradients, &
                                              displacement, result%hessian)
      result%has_hessian = .true.

      ! The undisplaced point supplies the energy and gradient that go with the
      ! Hessian.
      call method%calc_gradient(fragment, point)
      if (point%has_error) then
         result%error = point%error
         result%has_error = .true.
         call point%destroy()
         return
      end if
      result%energy = point%energy
      result%has_energy = point%has_energy
      if (point%has_gradient) then
         result%gradient = point%gradient
         result%has_gradient = .true.
      end if
      if (point%has_dipole) then
         result%dipole = point%dipole
         result%has_dipole = .true.
      end if
      call point%destroy()

      ! Dipole derivatives, hence IR intensities, only when every displacement
      ! returned a dipole; otherwise absent rather than silently zero.
      if (have_dipoles) then
         call finite_diff_dipole_derivatives(n_atoms, forward_dipoles, backward_dipoles, &
                                             displacement, result%dipole_derivatives)
         result%has_dipole_derivatives = .true.
         call report_translational_sum_rule(result%dipole_derivatives, fragment%charge, is_verbose)
      else
         result%has_dipole_derivatives = .false.
      end if
      deallocate (forward_dipoles, backward_dipoles)

      if (is_verbose) call logger%info("  Hessian assembled")
   end subroutine finite_difference_hessian

   subroutine finite_difference_hessian_energies(method, fragment, result, verbose, displacement_in)
      !! Hessian of `method` at `fragment`, by second differences of energies
      !!
      !! No gradient is used, so any method that can return an energy has a
      !! Hessian this way. `1 + 18 N**2` energy evaluations for `N` atoms (see
      !! `energy_hessian_task_count`), each a full SCF, so the energies need to
      !! be converged much tighter than for a gradient: the second difference
      !! divides their noise by `h**2`. The gradient comes from the same
      !! energies at no extra cost, and the dipole derivatives from the
      !! dipoles of the single displacements.
      class(qc_method_t), intent(in) :: method            !! Any method that returns an energy
      type(physical_fragment_t), intent(in) :: fragment   !! Reference geometry
      type(calculation_result_t), intent(out) :: result   !! Energy, gradient and Hessian
      logical, intent(in), optional :: verbose
      real(dp), intent(in), optional :: displacement_in
         !! Step in Bohr; `DEFAULT_DISPLACEMENT` when absent.

      type(physical_fragment_t) :: displaced
      type(calculation_result_t) :: point
      real(dp), allocatable :: energies(:)
      real(dp), allocatable :: forward_dipoles(:, :), backward_dipoles(:, :)
      real(dp), allocatable :: numerical_gradient(:, :)
      real(dp) :: displacement
      integer :: n_atoms, n_coords, n_tasks, task, coord_i, sign_i, coord_j, sign_j
      logical :: is_verbose, have_dipoles

      is_verbose = .false.
      if (present(verbose)) is_verbose = verbose

      n_atoms = fragment%n_atoms
      n_coords = 3*n_atoms
      n_tasks = energy_hessian_task_count(n_coords)
      displacement = DEFAULT_DISPLACEMENT
      if (present(displacement_in)) then
         if (displacement_in > 0.0_dp) displacement = displacement_in
      end if

      if (is_verbose) then
         call logger%info("Hessian by second differences of energies")
         call logger%info("  Atoms:                 "//to_char(n_atoms))
         call logger%info("  Energy evaluations:    "//to_char(n_tasks))
         call logger%info("  Displacement:          "//to_char(displacement)//" Bohr")
      end if

      allocate (energies(n_tasks))
      allocate (forward_dipoles(n_coords, 3), backward_dipoles(n_coords, 3))
      forward_dipoles = 0.0_dp
      backward_dipoles = 0.0_dp
      have_dipoles = .true.

      do task = 1, n_tasks
         call energy_hessian_task(task, n_coords, coord_i, sign_i, coord_j, sign_j)
         call energy_hessian_geometry(fragment, task, displacement, displaced)
         call method%calc_energy(displaced, point)
         call displaced%destroy()
         if (point%has_error .or. .not. point%has_energy) then
            call fail(result, point, displacement_kind(coord_i, coord_j), task, "energy")
            return
         end if
         energies(task) = point%energy%total()

         if (task == 1) then
            ! The undisplaced point supplies the energy that goes with the Hessian.
            result%energy = point%energy
            result%has_energy = .true.
            if (point%has_dipole) then
               result%dipole = point%dipole
               result%has_dipole = .true.
            end if
         else if (coord_j == 0) then
            ! Dipoles of the single displacements give dmu/dR; the paired ones add nothing.
            if (.not. point%has_dipole) then
               have_dipoles = .false.
            else if (sign_i > 0) then
               forward_dipoles(coord_i, :) = point%dipole
            else
               backward_dipoles(coord_i, :) = point%dipole
            end if
         end if
         call point%destroy()

         if (is_verbose .and. mod(task, max(1, n_tasks/10)) == 0) then
            call logger%info("  Completed "//to_char(task)//"/"//to_char(n_tasks)//" energy evaluations")
         end if
      end do

      call finite_diff_hessian_from_energies(n_coords, energies, displacement, result%hessian)
      result%has_hessian = .true.

      call finite_diff_gradient_from_energies(n_atoms, energies, displacement, numerical_gradient)
      call move_alloc(numerical_gradient, result%gradient)
      result%has_gradient = .true.

      ! Dipole derivatives, hence IR intensities, only when every single
      ! displacement returned a dipole; otherwise absent rather than zero.
      if (have_dipoles) then
         call finite_diff_dipole_derivatives(n_atoms, forward_dipoles, backward_dipoles, &
                                             displacement, result%dipole_derivatives)
         result%has_dipole_derivatives = .true.
         call report_translational_sum_rule(result%dipole_derivatives, fragment%charge, is_verbose)
      end if
      deallocate (forward_dipoles, backward_dipoles, energies)

      if (is_verbose) call logger%info("  Hessian assembled")
   end subroutine finite_difference_hessian_energies

   pure function displacement_kind(coord_i, coord_j) result(description)
      !! How to describe one energy evaluation of a Hessian from energies in an error
      integer, intent(in) :: coord_i, coord_j  !! As returned by `energy_hessian_task`
      character(len=:), allocatable :: description

      if (coord_i == 0) then
         description = "reference"
      else if (coord_j == 0) then
         description = "single"
      else
         description = "paired"
      end if
   end function displacement_kind

   subroutine suggest_numerical_hessian(result)
      !! Point a refusal at `keywords.hessian.method = "numerical"`
      !!
      !! A Hessian by differences of gradients needs the method's gradient, and
      !! several methods have none (coupled cluster, a continuum solvent on the
      !! CPU backend, a spin-scaled MP2). Their refusal is the gradient's, so it
      !! says to run an energy; this adds the way to a Hessian. Only a
      !! validation error is touched, since an SCF that failed to converge is
      !! not helped by a different route.
      type(calculation_result_t), intent(inout) :: result

      if (.not. result%has_error) return
      if (result%error%get_code() /= ERROR_VALIDATION) return
      call result%error%set(ERROR_VALIDATION, result%error%get_message()// &
                            " If it is the gradient that is not available, "// &
                            'keywords.hessian.method = "numerical" takes the Hessian '// &
                            "from energies alone and needs none.")
   end subroutine suggest_numerical_hessian

   subroutine report_translational_sum_rule(dipole_derivatives, charge, verbose)
      !! Check the translational sum rule on the dipole derivatives
      !!
      !!   sum_A  d(mu_k)/d(R_{A,k'})  =  Q * delta_{k,k'}
      !!
      !! A reference-free test of the dipole derivatives: a violation means
      !! dmu/dR is wrong, and shows up as spurious IR intensity on the
      !! translational modes.
      real(dp), intent(in) :: dipole_derivatives(:, :)  !! (3, 3*n_atoms)
      integer, intent(in) :: charge
      logical, intent(in) :: verbose

      real(dp) ::  residual
      real(dp) :: sums(3, 3)
      integer :: n_atoms, iatom, k, kp
      ! TODO(mqc): the warning branch below writes 100 characters into this
      ! 96-character buffer, so a sum rule that is actually violated raises an
      ! output-overflow runtime error instead of printing the warning.
      character(len=96) :: line

      n_atoms = size(dipole_derivatives, 2)/3
      sums = 0.0_dp
      do iatom = 1, n_atoms
         do kp = 1, 3
            do k = 1, 3
               sums(k, kp) = sums(k, kp) + dipole_derivatives(k, 3*(iatom - 1) + kp)
            end do
         end do
      end do

      ! Expected: charge on the diagonal, zero off it.
      residual = 0.0_dp
      do kp = 1, 3
         do k = 1, 3
            if (k == kp) then
               residual = max(residual, abs(sums(k, kp) - real(charge, dp)))
            else
               residual = max(residual, abs(sums(k, kp)))
            end if
         end do
      end do

      if (residual > 1.0e-4_dp) then
         write (line, "(a,es10.2,a)") &
            "IR: dipole-derivative translational sum rule off by ", residual, &
            " a.u. -- intensities may be unreliable"
         call logger%warning(trim(line))
      else if (verbose) then
         write (line, "(a,es10.2,a)") &
            "IR: dipole-derivative translational sum rule satisfied to ", residual, " a.u."
         call logger%info(trim(line))
      end if
   end subroutine report_translational_sum_rule

   subroutine fail(result, point, direction, index, quantity)
      !! Propagate a failed displacement outward with enough context to locate it
      type(calculation_result_t), intent(inout) :: result
      type(calculation_result_t), intent(inout) :: point
      character(len=*), intent(in) :: direction
      integer, intent(in) :: index
      character(len=*), intent(in), optional :: quantity  !! What was missing; "gradient" when absent

      character(len=:), allocatable :: missing

      missing = "gradient"
      if (present(quantity)) missing = quantity

      if (point%has_error) then
         result%error = point%error
         call result%error%add_context("Hessian: "//direction//" displacement "//to_char(index))
      else
         call result%error%set(ERROR_GENERIC, "Hessian: no "//missing//" returned for "// &
                               direction//" displacement "//to_char(index))
      end if
      result%has_error = .true.
      call point%destroy()
   end subroutine fail

end module mqc_semi_numerical_hessian
