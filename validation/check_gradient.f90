program check_gradient
   !! Compare the analytic cuEST gradient against central finite differences.
   !!
   !! This is the real test of the gradient: every term (Pulay, Hellmann-
   !! Feynman, density-fitted J/K, XC) has to be right and correctly signed
   !! for the two to agree, so a factor of two or a flipped sign shows up
   !! immediately. Both HF and DFT should match to roughly the SCF
   !! convergence: measured at 3.5e-8 (HF) and 3.4e-8 Ha/Bohr (PBE0) for
   !! water/def2-SVP, which is the finite-difference noise floor rather than
   !! any deficiency in the analytic expression.
   !!
   !! Usage:  ./check_gradient [functional]     (default: Hartree-Fock)
   use pic_types, only: dp
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_result_types, only: calculation_result_t
   use mqc_cuest_driver, only: cuest_scf_settings_t, run_cuest_scf
   implicit none

   type(cuest_scf_settings_t) :: settings
   type(physical_fragment_t) :: fragment
   type(calculation_result_t) :: result
   real(dp), allocatable :: analytic(:, :), numeric(:, :)
   real(dp) :: e_plus, e_minus, step, max_dev, rms
   integer :: iatom, ixyz, n_atoms, n_arg
   character(len=32) :: functional

   functional = ''
   n_arg = command_argument_count()
   if (n_arg >= 1) call get_command_argument(1, functional)

   step = 5.0e-4_dp   ! Bohr; small enough for truncation, large enough vs SCF noise

   call build_water(fragment)
   n_atoms = fragment%n_atoms

   settings%basis_set = 'def2-svp'
   settings%aux_basis_set = 'def2-universal-jkfit'
   settings%functional = functional
   ! The gradient inherits roughly the square root of the density error, and
   ! the finite differences need energies far tighter than the step size, so
   ! both tolerances are pushed well past the defaults.
   settings%energy_tol = 1.0e-11_dp
   settings%density_tol = 1.0e-8_dp
   settings%max_iter = 200

   if (len_trim(functional) == 0) then
      write (*, '(A)') "Gradient check: Hartree-Fock / def2-SVP / water"
   else
      write (*, '(A,A,A)') "Gradient check: ", trim(functional), " / def2-SVP / water"
   end if
   write (*, '(A,ES9.2,A)') "  finite-difference step: ", step, " Bohr (central)"

   ! ---- analytic ----------------------------------------------------------
   call run_cuest_scf(settings, fragment, result, want_gradient=.true.)
   if (result%has_error) then
      write (*, '(A)') "FAILED: "//result%error%get_message()
      error stop 1
   end if
   allocate (analytic(3, n_atoms))
   analytic = result%gradient
   write (*, '(A,F20.12)') "  energy: ", result%energy%scf

   ! ---- central differences ------------------------------------------------
   allocate (numeric(3, n_atoms))
   do iatom = 1, n_atoms
      do ixyz = 1, 3
         fragment%coordinates(ixyz, iatom) = fragment%coordinates(ixyz, iatom) + step
         e_plus = single_point(settings, fragment)
         fragment%coordinates(ixyz, iatom) = fragment%coordinates(ixyz, iatom) - 2.0_dp*step
         e_minus = single_point(settings, fragment)
         fragment%coordinates(ixyz, iatom) = fragment%coordinates(ixyz, iatom) + step
         numeric(ixyz, iatom) = (e_plus - e_minus)/(2.0_dp*step)
      end do
   end do

   ! ---- report --------------------------------------------------------------
   write (*, '(A)') ""
   write (*, '(A)') "  atom  cmp        analytic          finite-diff         difference"
   do iatom = 1, n_atoms
      do ixyz = 1, 3
         write (*, '(2X,I4,3X,I2,3X,3F20.12)') iatom, ixyz, &
            analytic(ixyz, iatom), numeric(ixyz, iatom), &
            analytic(ixyz, iatom) - numeric(ixyz, iatom)
      end do
   end do

   max_dev = maxval(abs(analytic - numeric))
   rms = sqrt(sum((analytic - numeric)**2)/real(3*n_atoms, dp))
   write (*, '(A)') ""
   write (*, '(A,ES12.4,A)') "  max |analytic - numeric| : ", max_dev, " Ha/Bohr"
   write (*, '(A,ES12.4,A)') "  rms deviation            : ", rms, " Ha/Bohr"
   write (*, '(A,ES12.4,A)') "  gradient norm            : ", sqrt(sum(analytic**2)), " Ha/Bohr"

contains

   function single_point(settings, fragment) result(energy)
      type(cuest_scf_settings_t), intent(in) :: settings
      type(physical_fragment_t), intent(in) :: fragment
      real(dp) :: energy
      type(calculation_result_t) :: point

      call run_cuest_scf(settings, fragment, point)
      if (point%has_error) then
         write (*, '(A)') "displaced point FAILED: "//point%error%get_message()
         error stop 1
      end if
      energy = point%energy%scf
   end function single_point

   subroutine build_water(fragment)
      !! Water in Bohr, deliberately off equilibrium so the gradient is large
      !! enough that a sign error or a missing term cannot hide in the noise.
      type(physical_fragment_t), intent(out) :: fragment

      fragment%n_atoms = 3
      allocate (fragment%element_numbers(3), fragment%coordinates(3, 3))
      fragment%element_numbers = [8, 1, 1]
      fragment%coordinates(:, 1) = [0.0_dp, 0.0_dp, 0.19043_dp]
      fragment%coordinates(:, 2) = [0.0_dp, 1.51983_dp, -0.88401_dp]
      fragment%coordinates(:, 3) = [0.0_dp, -1.41983_dp, -0.94401_dp]
      fragment%charge = 0
      fragment%multiplicity = 1
      fragment%nelec = 10
      fragment%n_caps = 0
   end subroutine build_water

end program check_gradient
