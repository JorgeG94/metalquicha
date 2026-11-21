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
   use mqc_czt_gradient, only: czt_scf_gradient
   use mqc_czt_multipole, only: multipole_matrices, dipole_integral_derivatives
   use mqc_czt_xc, only: xc_context_t, xc_context_create
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_scf_convergence, only: scf_convergence_t
   implicit none
   private

   public :: polarizability_derivatives
   public :: polarizability_derivatives_field

   ! The displacement is the Hessian's default, 5e-3 Bohr. A central
   ! difference's error is quadratic in it: 3e-5 Bohr^2 on water's
   ! derivatives of order ten, measured by the translational sum rule, which
   ! is the sixth figure of an activity.
   real(dp), parameter :: RAMAN_FIELD = 1.0e-3_dp
      !! Field strength in atomic units for the finite-field route. The error
      !! of the second difference is quadratic in it and was measured on
      !! water against the displacement route: 2e-3 Bohr^2 at 8e-3, 5e-4 at
      !! 4e-3, 1e-4 at 2e-3, 3e-5 here, a clean factor of four each time. At
      !! half of this the gradient's own noise, divided by the field squared,
      !! begins to show.
   integer, parameter :: N_FIELDS = 12
      !! Fields the second differences need: plus and minus each axis, and
      !! plus and minus each pair of axes together

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

   subroutine polarizability_derivatives_field(settings, mol, nelec, reference_density, &
                                               accelerator, diis_vectors, convergence, &
                                               alpha, dalpha, error, xc)
      !! `d alpha / d x` from the analytic gradient in a finite electric field
      !!
      !! In a uniform field `F` the energy is `E0 - mu.F - (1/2) F alpha F`,
      !! so the polarizability's nuclear derivative is minus the second field
      !! derivative of the nuclear gradient `g`:
      !!
      !!     d alpha_kl / d x = - d2 g_x / dF_k dF_l
      !!
      !! taken by second differences over twelve fields and the field-free
      !! gradient. **Thirteen self-consistent fields and gradients whatever
      !! the size of the molecule**, where the displacement route needs six
      !! per atom, and every one of them at the reference geometry, so the
      !! molecule and its quadrature grid are built once.
      !!
      !! The field enters the SCF as the one-electron operator `+ r.F`. The
      !! gradient of that term is added here, since the SCF gradient does not
      !! know about it: the basis functions' part,
      !! `sum_a F_a Tr(D d r_a / d R)`, and the nuclei's, `- Z F`.
      type(cuest_scf_settings_t), intent(in) :: settings
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: nelec
      real(dp), intent(in) :: reference_density(:, :)   !! Total, field free
      integer, intent(in) :: accelerator, diis_vectors
      type(scf_convergence_t), intent(in) :: convergence
      real(dp), intent(out) :: alpha(3, 3)   !! Field free, Bohr^3
      real(dp), allocatable, intent(out) :: dalpha(:, :)
         !! (6, 3*n_atoms) in Bohr^2, laid out as `polarizability_derivatives`
         !! lays it out
      type(error_t), intent(inout) :: error
      type(xc_context_t), intent(inout), optional :: xc
         !! The reference's functional. Absent or inactive is Hartree-Fock.

      integer, parameter :: PAIRS(2, 3) = reshape([1, 2, 1, 3, 2, 3], [2, 3])
      real(dp), allocatable :: dip(:, :, :), ddip(:, :, :, :, :)
      real(dp), allocatable :: g0(:, :), g(:, :, :), second(:, :, :)
      real(dp) :: origin(3), fields(3, N_FIELDS)
      type(scf_convergence_t) :: tight
      type(rhf_result_t) :: scf
      integer :: natm, k, l, ipair, ifield, iatom, comp, col
      logical :: kohn_sham

      alpha = 0.0_dp
      if (error%has_error()) return
      natm = mol%natm
      kohn_sham = .false.
      if (present(xc)) kohn_sham = xc%active
      tight = convergence%for_differencing()

      ! The origin `r` is measured from. It moves the energy of a charged
      ! molecule by a constant per unit field and its gradient by a term
      ! linear in the field, neither of which a second difference keeps.
      origin = sum(mol%coords, 2)/real(natm, dp)
      call multipole_matrices(mol, origin, 1, dip, error)
      if (error%has_error()) return
      call dipole_integral_derivatives(mol, origin, ddip, error)
      if (error%has_error()) return

      ! Fields 1-6: plus and minus each axis. Fields 7-12: plus and minus
      ! each pair of axes together.
      fields = 0.0_dp
      do k = 1, 3
         fields(k, 2*k - 1) = RAMAN_FIELD
         fields(k, 2*k) = -RAMAN_FIELD
      end do
      do ipair = 1, 3
         fields(PAIRS(:, ipair), 5 + 2*ipair) = RAMAN_FIELD
         fields(PAIRS(:, ipair), 6 + 2*ipair) = -RAMAN_FIELD
      end do

      call logger%info("  Raman: differencing the gradient over "//to_char(N_FIELDS)// &
                       " electric fields")
      call field_gradient([0.0_dp, 0.0_dp, 0.0_dp], g0, keep=scf)
      if (error%has_error()) return
      if (kohn_sham) then
         call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                    scf%n_occupied, alpha, error, tol=PROPERTY_RESPONSE_TOL, &
                                    xc=xc, density=scf%density)
      else
         call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                    scf%n_occupied, alpha, error, tol=PROPERTY_RESPONSE_TOL)
      end if
      if (error%has_error()) return

      allocate (g(3, natm, N_FIELDS))
      do ifield = 1, N_FIELDS
         block
            real(dp), allocatable :: one(:, :)

            call field_gradient(fields(:, ifield), one)
            if (error%has_error()) then
               call error%add_context("Raman: a gradient in a finite field")
               return
            end if
            g(:, :, ifield) = one
         end block
      end do

      ! Second derivatives of the gradient: the diagonal ones from each
      ! axis, the mixed ones from the pair fields less the two diagonals.
      allocate (second(3, natm, N_SYMMETRIC))
      do k = 1, 3
         second(:, :, k) = (g(:, :, 2*k - 1) + g(:, :, 2*k) - 2.0_dp*g0)/RAMAN_FIELD**2
      end do
      do ipair = 1, 3
         k = PAIRS(1, ipair)
         l = PAIRS(2, ipair)
         second(:, :, 3 + ipair) = 0.5_dp*((g(:, :, 5 + 2*ipair) + g(:, :, 6 + 2*ipair) &
                                            - 2.0_dp*g0)/RAMAN_FIELD**2 &
                                           - second(:, :, k) - second(:, :, l))
      end do

      ! Components xx, yy, zz, xy, xz, yz, which is the order of `PAIRS`.
      allocate (dalpha(N_SYMMETRIC, 3*natm))
      do iatom = 1, natm
         do comp = 1, 3
            col = 3*(iatom - 1) + comp
            dalpha(:, col) = -second(comp, iatom, :)
         end do
      end do

   contains

      subroutine field_gradient(field, gradient, keep)
         !! The nuclear gradient with the molecule in a uniform field
         real(dp), intent(in) :: field(3)
         real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, natm)
         type(rhf_result_t), intent(out), optional :: keep
            !! The SCF it was taken from

         type(rhf_result_t) :: here
         real(dp), allocatable :: h_field(:, :)
         integer :: a, b, x

         allocate (h_field(mol%nao, mol%nao))
         h_field = 0.0_dp
         do a = 1, 3
            h_field = h_field + field(a)*dip(:, :, a)
         end do
         call run_czt_rhf(mol, nelec, settings%max_iter, DIFFERENCED_SCF_ENERGY_TOL, &
                          DIFFERENCED_SCF_DENSITY_TOL, .false., here, error, &
                          guess=SCF_GUESS_SAD, &
                          guess_density=reference_density, xc=xc, h_extra=h_field, &
                          level_shift=settings%level_shift, &
                          linear_dependence=settings%linear_dependence, &
                          accelerator=accelerator, diis_vectors=diis_vectors, &
                          incremental_fock=settings%incremental_fock, &
                          grad_tol=tight%gradient_tolerance, convergence=tight)
         if (error%has_error()) return
         if (.not. here%converged) then
            call error%set(ERROR_VALIDATION, "the reference did not converge in a "// &
                           "finite field")
            return
         end if
         if (kohn_sham) then
            call czt_scf_gradient(mol, here%density, orbitals=here%orbitals, &
                                  orbital_energies=here%orbital_energies, &
                                  n_occupied=here%n_occupied, gradient=gradient, &
                                  error=error, xc=xc)
         else
            call czt_scf_gradient(mol, here%density, orbitals=here%orbitals, &
                                  orbital_energies=here%orbital_energies, &
                                  n_occupied=here%n_occupied, gradient=gradient, &
                                  error=error)
         end if
         if (error%has_error()) return
         ! The field term's own derivative: the basis functions carry `r`
         ! with them, and each nucleus feels `-Z F` directly.
         do x = 1, mol%natm
            do b = 1, 3
               do a = 1, 3
                  gradient(b, x) = gradient(b, x) &
                                   + field(a)*sum(here%density*ddip(:, :, a, b, x))
               end do
               gradient(b, x) = gradient(b, x) - mol%charges(x)*field(b)
            end do
         end do
         if (present(keep)) keep = here
      end subroutine field_gradient
   end subroutine polarizability_derivatives_field

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
