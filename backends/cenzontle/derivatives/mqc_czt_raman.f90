!! Cartesian derivatives of the static polarizability, for Raman intensities
module mqc_czt_raman
   !! A Raman line is as strong as the polarizability changes along its
   !! normal mode. This module supplies the Cartesian half of that: the
   !! derivative of the static polarizability tensor with respect to each
   !! nuclear coordinate, which `mqc_vibrational_analysis` then projects onto
   !! the modes.
   !!
   !! **By central differences of an analytic polarizability.** Each
   !! coordinate is displaced both ways, the closed-shell reference is
   !! reconverged, and the coupled-perturbed polarizability is taken there:
   !! `6 N` self-consistent fields and as many response solves. The analytic
   !! route is a third derivative of the energy and is not built.
   !!
   !! **Static, so off resonance.** The polarizability is the zero-frequency
   !! one. A measured Raman spectrum is taken at a laser frequency, and near
   !! an electronic absorption the two differ a great deal.
   use pic_types, only: dp
   use pic_io, only: to_char
   use pic_logger, only: logger => global_logger
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, PROPERTY_RESPONSE_TOL, &
                                       DIFFERENCED_SCF_ENERGY_TOL, &
                                       DIFFERENCED_SCF_DENSITY_TOL
   use mqc_symmetric_tensor, only: pack_symmetric, N_SYMMETRIC
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, SCF_GUESS_SAD
   use mqc_czt_cphf, only: static_polarizability
   use mqc_czt_xc, only: xc_context_t, xc_context_create
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_scf_convergence, only: scf_convergence_t
   implicit none
   private

   public :: polarizability_derivatives

   ! The displacement is the Hessian's default, 5e-3 Bohr. A central
   ! difference's error is quadratic in it: 3e-5 Bohr^2 on water's
   ! derivatives of order ten, measured by the translational sum rule, which
   ! is the sixth figure of an activity.

contains

   subroutine polarizability_derivatives(settings, element_numbers, symbols, coords, &
                                         nelec, reference_density, accelerator, &
                                         diis_vectors, convergence, alpha, dalpha, error)
      !! `d alpha / d x` for every Cartesian coordinate, and `alpha` itself
      !!
      !! Closed-shell Hartree-Fock or Kohn-Sham, exact integrals. The
      !! reference density seeds every displaced SCF.
      type(cuest_scf_settings_t), intent(in) :: settings
         !! Read: the basis and its Cartesian flag, the functional and its
         !! grid, and the SCF iteration limit
      integer, intent(in) :: element_numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)   !! (3, n_atoms), Bohr
      integer, intent(in) :: nelec
      real(dp), intent(in) :: reference_density(:, :)   !! Total, at `coords`
      integer, intent(in) :: accelerator, diis_vectors
         !! As the reference SCF resolved them from the deck
      type(scf_convergence_t), intent(in) :: convergence
         !! The reference's rule. Its metric is kept and its tolerances are
         !! tightened to what a differenced polarizability needs.
      real(dp), intent(out) :: alpha(3, 3)   !! At `coords`, Bohr^3
      real(dp), allocatable, intent(out) :: dalpha(:, :)
         !! (N_SYMMETRIC, 3*n_atoms) in Bohr^2: the components
         !! `xx, yy, zz, xy, xz, yz` of the derivative with respect to
         !! coordinate `3*(atom-1) + xyz`
      type(error_t), intent(inout) :: error

      real(dp) :: moved(3, size(coords, 2)), plus(3, 3), minus(3, 3), diff(3, 3)
      type(scf_convergence_t) :: tight
      integer :: n_atoms, iatom, comp

      if (error%has_error()) return
      n_atoms = size(coords, 2)
      tight = convergence%for_differencing()
      allocate (dalpha(N_SYMMETRIC, 3*n_atoms))
      dalpha = 0.0_dp

      call polarizability_at(settings, element_numbers, symbols, coords, nelec, &
                             reference_density, accelerator, diis_vectors, tight, &
                             alpha, error)
      if (error%has_error()) return

      call logger%info("  Raman: differencing the polarizability over "// &
                       to_char(2*3*n_atoms)//" displaced geometries")
      do iatom = 1, n_atoms
         do comp = 1, 3
            moved = coords
            moved(comp, iatom) = coords(comp, iatom) + DEFAULT_DISPLACEMENT
            call polarizability_at(settings, element_numbers, symbols, moved, nelec, &
                                   reference_density, accelerator, diis_vectors, tight, &
                                   plus, error)
            moved(comp, iatom) = coords(comp, iatom) - DEFAULT_DISPLACEMENT
            call polarizability_at(settings, element_numbers, symbols, moved, nelec, &
                                   reference_density, accelerator, diis_vectors, tight, &
                                   minus, error)
            if (error%has_error()) then
               call error%add_context("Raman: a displaced polarizability")
               return
            end if
            diff = (plus - minus)/(2.0_dp*DEFAULT_DISPLACEMENT)
            dalpha(:, 3*(iatom - 1) + comp) = pack_symmetric(diff)
         end do
      end do
   end subroutine polarizability_derivatives

   subroutine polarizability_at(settings, element_numbers, symbols, coords, nelec, &
                                guess_density, accelerator, diis_vectors, convergence, &
                                alpha, error)
      !! Converge the reference at one geometry and return its polarizability
      type(cuest_scf_settings_t), intent(in) :: settings
      integer, intent(in) :: element_numbers(:)
      character(len=*), intent(in) :: symbols(:)
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: nelec
      real(dp), intent(in) :: guess_density(:, :)
      integer, intent(in) :: accelerator, diis_vectors
      type(scf_convergence_t), intent(in) :: convergence
      real(dp), intent(out) :: alpha(3, 3)
      type(error_t), intent(inout) :: error

      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: xc
      logical :: kohn_sham

      alpha = 0.0_dp
      if (error%has_error()) return
      call build_czt_molecule(element_numbers, symbols, coords, trim(settings%basis_set), &
                              mol, error, force_cartesian=settings%cartesian)
      if (error%has_error()) return

      ! An inactive `xc` is Hartree-Fock to the SCF and the response alike.
      kohn_sham = len_trim(settings%functional) > 0
      if (kohn_sham) then
         call xc_context_create(mol, trim(settings%functional), xc, error, &
                                level=settings%grid_level, polarized=.false., &
                                nlc_level=settings%nlc_grid_level, &
                                screen_tol=settings%screening_tolerance, &
                                point_block=settings%block_size, &
                                small_rho_cutoff=settings%small_rho_cutoff, &
                                n_radial=settings%radial_points, &
                                n_angular=settings%angular_points)
         if (error%has_error()) then
            call mol%destroy()
            return
         end if
      end if
      call run_czt_rhf(mol, nelec, settings%max_iter, DIFFERENCED_SCF_ENERGY_TOL, &
                       DIFFERENCED_SCF_DENSITY_TOL, .false., scf, error, &
                       guess=SCF_GUESS_SAD, guess_density=guess_density, xc=xc, &
                       level_shift=settings%level_shift, &
                       linear_dependence=settings%linear_dependence, &
                       accelerator=accelerator, diis_vectors=diis_vectors, &
                       incremental_fock=settings%incremental_fock, &
                       grad_tol=convergence%gradient_tolerance, &
                       convergence=convergence)
      if (.not. error%has_error() .and. .not. scf%converged) then
         call error%set(ERROR_VALIDATION, "the reference did not converge at a "// &
                        "displaced geometry")
      end if
      if (.not. error%has_error()) then
         if (kohn_sham) then
            call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                       scf%n_occupied, alpha, error, &
                                       tol=PROPERTY_RESPONSE_TOL, xc=xc, &
                                       density=scf%density)
         else
            call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                       scf%n_occupied, alpha, error, &
                                       tol=PROPERTY_RESPONSE_TOL)
         end if
      end if
      if (kohn_sham) call xc%destroy()
      call mol%destroy()
   end subroutine polarizability_at

end module mqc_czt_raman
