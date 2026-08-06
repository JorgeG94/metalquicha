!! Hessians by finite differences of analytic gradients
module mqc_semi_numerical_hessian
   !! Builds a Hessian from central differences of a method's analytic
   !! gradients:
   !!
   !!   H[i,j] = (g_j(x_i + h) - g_j(x_i - h)) / 2h
   !!
   !! This is *semi*-numerical: only one derivative is taken numerically, the
   !! other is analytic. That matters for accuracy -- differencing energies
   !! twice would lose roughly half the available digits, whereas differencing
   !! gradients once keeps most of them.
   !!
   !! The routine is written against `qc_method_t` rather than any concrete
   !! method, so anything with a working `calc_gradient` gets a Hessian for
   !! free.
   !!
   !! Cost is 6N gradient evaluations for N atoms, and each of those is a full
   !! SCF. The displaced geometries are independent, so this is the obvious
   !! place to distribute work across ranks -- metalquicha already does that
   !! for the unfragmented Hessian path.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: ERROR_GENERIC
   use mqc_finite_differences, only: generate_perturbed_geometries, displaced_geometry_t, &
                                     finite_diff_hessian_from_gradients, &
                                     finite_diff_dipole_derivatives, DEFAULT_DISPLACEMENT
   implicit none
   private

   public :: finite_difference_hessian  !! Hessian from central differences of gradients

contains

   subroutine finite_difference_hessian(method, fragment, result, verbose)
      !! Hessian of `method` at `fragment`, by central differences of gradients
      class(qc_method_t), intent(in) :: method            !! Any method with analytic gradients
      type(physical_fragment_t), intent(in) :: fragment   !! Reference geometry
      type(calculation_result_t), intent(out) :: result   !! Energy, gradient and Hessian
      logical, intent(in), optional :: verbose

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

      if (is_verbose) then
         call logger%info("Hessian by central differences of analytic gradients")
         call logger%info("  Atoms:                 "//to_char(n_atoms))
         call logger%info("  Gradient evaluations:  "//to_char(2*n_displacements))
         call logger%info("  Displacement:          "//to_char(displacement)//" Bohr")
      end if

      call generate_perturbed_geometries(fragment, displacement, forward_geoms, backward_geoms)
      allocate (forward_gradients(n_displacements, 3, n_atoms))
      allocate (backward_gradients(n_displacements, 3, n_atoms))
      ! The same displacements give dmu/dR, and so IR intensities, for nothing
      ! beyond storing the dipole each method already returns.
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
      ! Hessian. It is one extra SCF on top of the 6N, and skipping it would
      ! mean reporting a Hessian with no energy beside it.
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

      ! Dipole derivatives, hence IR intensities, whenever every displacement
      ! returned a dipole. A method that does not supply one leaves them
      ! absent rather than silently zero.
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

   subroutine report_translational_sum_rule(dipole_derivatives, charge, verbose)
      !! Check the translational sum rule on the dipole derivatives
      !!
      !!   sum_A  d(mu_k)/d(R_{A,k'})  =  Q * delta_{k,k'}
      !!
      !! Rigid translation of the whole molecule shifts every nucleus and the
      !! whole electron cloud together, so the dipole changes only by Q times
      !! the shift -- zero for a neutral system. It is a reference-free test of
      !! the dipole derivatives: a violation means dmu/dR is wrong, and shows up
      !! as spurious IR intensity on the translational modes.
      real(dp), intent(in) :: dipole_derivatives(:, :)  !! (3, 3*n_atoms)
      integer, intent(in) :: charge
      logical, intent(in) :: verbose

      real(dp) ::  residual
      real(dp) :: sums(3, 3)
      integer :: n_atoms, iatom, k, kp
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

      ! Expected: charge on the diagonal, zero off it. Total nuclear charge and
      ! electron count balance to `charge` (in atomic units the dipole gradient
      ! sum rule uses the total molecular charge).
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

   subroutine fail(result, point, direction, index)
      !! Propagate a failed displacement outward with enough context to locate it
      type(calculation_result_t), intent(inout) :: result
      type(calculation_result_t), intent(inout) :: point
      character(len=*), intent(in) :: direction
      integer, intent(in) :: index

      if (point%has_error) then
         result%error = point%error
         call result%error%add_context("Hessian: "//direction//" displacement "//to_char(index))
      else
         call result%error%set(ERROR_GENERIC, "Hessian: no gradient returned for "// &
                               direction//" displacement "//to_char(index))
      end if
      result%has_error = .true.
      call point%destroy()
   end subroutine fail

end module mqc_semi_numerical_hessian
