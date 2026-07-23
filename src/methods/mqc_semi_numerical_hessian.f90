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
                                     finite_diff_hessian_from_gradients, DEFAULT_DISPLACEMENT
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
      type(calculation_result_t) :: point
      real(dp) :: displacement
      integer :: n_atoms, n_displacements, i
      logical :: is_verbose

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

      do i = 1, n_displacements
         call method%calc_gradient(forward_geoms(i)%geometry, point)
         if (point%has_error .or. .not. point%has_gradient) then
            call fail(result, point, "forward", i)
            return
         end if
         forward_gradients(i, :, :) = point%gradient
         call point%destroy()

         call method%calc_gradient(backward_geoms(i)%geometry, point)
         if (point%has_error .or. .not. point%has_gradient) then
            call fail(result, point, "backward", i)
            return
         end if
         backward_gradients(i, :, :) = point%gradient
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
      call point%destroy()

      ! No dipole derivatives: the cuEST backend does not compute dipoles yet,
      ! so IR intensities are unavailable rather than silently zero.
      result%has_dipole_derivatives = .false.

      if (is_verbose) call logger%info("  Hessian assembled")
   end subroutine finite_difference_hessian

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
