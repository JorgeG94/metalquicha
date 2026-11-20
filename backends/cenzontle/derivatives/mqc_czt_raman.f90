!! Cartesian derivatives of the static polarizability, for Raman intensities
module mqc_czt_raman
   !! A Raman line is as strong as the polarizability changes along its
   !! normal mode. This module supplies the Cartesian half of that: the
   !! derivative of the static polarizability tensor with respect to each
   !! nuclear coordinate, which `mqc_vibrational_analysis` then projects onto
   !! the modes.
   !!
   !! **Three routes to one third derivative of the energy.**
   !! `polarizability_derivatives_analytic` differentiates the gradient twice
   !! in the field in closed form, for Hartree-Fock.
   !! `polarizability_derivatives_field` takes the same two derivatives by
   !! differences over twelve fields, and `polarizability_derivatives`
   !! differences the polarizability over `6 N` displaced geometries; both
   !! run under a functional too.
   !!
   !! **Static, so off resonance.** The polarizability is the zero-frequency
   !! one. A measured Raman spectrum is taken at a laser frequency, and near
   !! an electronic absorption the two differ a great deal.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, SCF_GUESS_SAD
   use pic_blas_interfaces, only: pic_gemm
   use mqc_czt_integrals, only: atom_ao_blocks
   use mqc_czt_cphf, only: static_polarizability, cphf_solve
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_response_product, only: response_mean_field
   use mqc_czt_gradient, only: czt_scf_gradient, one_electron_deriv, &
                               two_electron_gradient, hellmann_feynman_gradient, &
                               DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_czt_multipole, only: multipole_matrices, dipole_integral_derivatives
   use mqc_czt_xc, only: xc_context_t, xc_context_create
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_scf_convergence, only: scf_convergence_t
   implicit none
   private

   public :: polarizability_derivatives
   public :: polarizability_derivatives_field
   public :: polarizability_derivatives_analytic
   public :: RAMAN_FIELD
   public :: RAMAN_STEP

   real(dp), parameter :: RAMAN_STEP = 5.0e-3_dp
      !! Displacement in Bohr. The error of a central difference is quadratic
      !! in it: 3e-5 Bohr^2 on water's derivatives of order ten, measured by
      !! the translational sum rule, which is the sixth figure of an activity.
   real(dp), parameter :: SCF_ENERGY_TOL = 1.0e-12_dp
   real(dp), parameter :: SCF_DENSITY_TOL = 1.0e-10_dp
   real(dp), parameter :: SCF_GRAD_TOL = 1.0e-10_dp
      !! Far tighter than an energy needs: the polarizabilities are
      !! differenced, which divides their noise by the step. An orbital
      !! gradient of 1e-8 left the derivatives 3e-6 off.
   real(dp), parameter :: RESPONSE_TOL = 1.0e-10_dp
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
   integer, parameter :: N_COMPONENTS = 6
      !! Independent components of the symmetric tensor: xx, yy, zz, xy, xz, yz
   integer, parameter :: LINE_LEN = 160

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
         !! (6, 3*n_atoms) in Bohr^2: the components `xx, yy, zz, xy, xz, yz`
         !! of the derivative with respect to coordinate `3*(atom-1) + xyz`
      type(error_t), intent(inout) :: error

      real(dp) :: moved(3, size(coords, 2)), plus(3, 3), minus(3, 3), diff(3, 3)
      type(scf_convergence_t) :: tight
      integer :: n_atoms, iatom, comp, k
      character(len=LINE_LEN) :: line

      if (error%has_error()) return
      n_atoms = size(coords, 2)
      tight = convergence
      tight%tolerance = min(convergence%tolerance, SCF_ENERGY_TOL)
      tight%gradient_tolerance = SCF_GRAD_TOL
      allocate (dalpha(6, 3*n_atoms))
      dalpha = 0.0_dp

      call polarizability_at(settings, element_numbers, symbols, coords, nelec, &
                             reference_density, accelerator, diis_vectors, tight, &
                             alpha, error)
      if (error%has_error()) return

      write (line, "(a,i0,a)") "  Raman: differencing the polarizability over ", &
         6*n_atoms, " displaced geometries"
      call logger%info(trim(line))
      do iatom = 1, n_atoms
         do comp = 1, 3
            moved = coords
            moved(comp, iatom) = coords(comp, iatom) + RAMAN_STEP
            call polarizability_at(settings, element_numbers, symbols, moved, nelec, &
                                   reference_density, accelerator, diis_vectors, tight, &
                                   plus, error)
            moved(comp, iatom) = coords(comp, iatom) - RAMAN_STEP
            call polarizability_at(settings, element_numbers, symbols, moved, nelec, &
                                   reference_density, accelerator, diis_vectors, tight, &
                                   minus, error)
            if (error%has_error()) then
               call error%add_context("Raman: a displaced polarizability")
               return
            end if
            diff = (plus - minus)/(2.0_dp*RAMAN_STEP)
            ! The tensor is symmetric; the average of the two off-diagonal
            ! elements removes what the response solver's tolerance leaves.
            k = 3*(iatom - 1) + comp
            dalpha(1, k) = diff(1, 1)
            dalpha(2, k) = diff(2, 2)
            dalpha(3, k) = diff(3, 3)
            dalpha(4, k) = 0.5_dp*(diff(1, 2) + diff(2, 1))
            dalpha(5, k) = 0.5_dp*(diff(1, 3) + diff(3, 1))
            dalpha(6, k) = 0.5_dp*(diff(2, 3) + diff(3, 2))
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
      character(len=LINE_LEN) :: line

      alpha = 0.0_dp
      if (error%has_error()) return
      natm = mol%natm
      kohn_sham = .false.
      if (present(xc)) kohn_sham = xc%active
      tight = convergence
      tight%tolerance = min(convergence%tolerance, SCF_ENERGY_TOL)
      tight%gradient_tolerance = SCF_GRAD_TOL

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

      write (line, "(a,i0,a)") "  Raman: differencing the gradient over ", N_FIELDS, &
         " electric fields"
      call logger%info(trim(line))
      call field_gradient([0.0_dp, 0.0_dp, 0.0_dp], g0, keep=scf)
      if (error%has_error()) return
      if (kohn_sham) then
         call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                    scf%n_occupied, alpha, error, tol=RESPONSE_TOL, &
                                    xc=xc, density=scf%density)
      else
         call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                    scf%n_occupied, alpha, error, tol=RESPONSE_TOL)
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
      allocate (second(3, natm, 6))
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
      allocate (dalpha(6, 3*natm))
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
         call run_czt_rhf(mol, nelec, settings%max_iter, SCF_ENERGY_TOL, SCF_DENSITY_TOL, &
                          .false., here, error, guess=SCF_GUESS_SAD, &
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

   subroutine polarizability_derivatives_analytic(mol, orbitals, orbital_energies, n_occ, &
                                                  alpha, dalpha, error)
      !! `d alpha / d x` with no finite difference, for closed-shell Hartree-Fock
      !!
      !! The nuclear gradient holds in a uniform field `F` as it does without
      !! one,
      !!
      !!     g_x = Tr(D h^x) + (1/2) Tr(D G^x[D]) - Tr(W S^x) + F . Tr(D r^x)
      !!
      !! and the polarizability's derivative is minus its second field
      !! derivative. The integrals do not depend on the field, so that is the
      !! same contraction at the density's and the energy-weighted density's
      !! field derivatives:
      !!
      !!     d alpha_kl / d x = - [ Tr(D^kl h^x) - Tr(W^kl S^x)
      !!                            + Tr(D^kl G^x[D]) + Tr(D^k G^x[D^l])
      !!                            + Tr(D^k r_l^x) + Tr(D^l r_k^x) ]
      !!
      !! `D^k` comes from the three coupled-perturbed solves the
      !! polarizability itself needs and `D^kl` from six more, against the
      !! same operator with a right-hand side built from the first-order
      !! ones. No nuclear response is solved, so the cost does not grow with
      !! the number of atoms beyond that of a gradient: nine response solves
      !! and eighteen two-electron gradient contractions.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)          !! Canonical, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: alpha(3, 3)   !! Bohr^3
      real(dp), allocatable, intent(out) :: dalpha(:, :)
         !! (6, 3*n_atoms) in Bohr^2, laid out as `polarizability_derivatives`
         !! lays it out
      type(error_t), intent(inout) :: error

      integer, parameter :: PAIRS(2, N_COMPONENTS) = &
                            reshape([1, 1, 2, 2, 3, 3, 1, 2, 1, 3, 2, 3], [2, N_COMPONENTS])
      real(dp), allocatable :: dip(:, :, :), ddip(:, :, :, :, :)
      real(dp), allocatable :: u1(:, :, :), u2(:, :, :), rhs(:, :, :)
      real(dp), allocatable :: d0(:, :), d1(:, :, :), d2(:, :, :), w2(:, :, :)
      real(dp), allocatable :: dm1(:, :, :), dm2(:, :, :), f1(:, :, :), f2(:, :, :)
      real(dp), allocatable :: dm0(:, :), f0(:, :), wm(:, :)
      real(dp), allocatable :: g(:, :, :), scratch(:, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :)
      real(dp), allocatable :: q1(:, :, :), qpair(:, :), plus(:, :), minus(:, :)
      real(dp), allocatable :: second(:, :), trial(:, :)
      integer, allocatable :: offsets(:), counts(:)
      real(dp) :: origin(3)
      integer :: n_ao, n_mo, n_vir, natm, k, l, c, i, iatom, comp, p0, p1

      alpha = 0.0_dp
      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ
      natm = mol%natm

      origin = sum(mol%coords, 2)/real(natm, dp)
      call multipole_matrices(mol, origin, 1, dip, error)
      if (error%has_error()) return
      call dipole_integral_derivatives(mol, origin, ddip, error)
      if (error%has_error()) return
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      allocate (zero_h(n_ao, n_ao))
      zero_h = 0.0_dp

      call logger%info("  Raman: analytic polarizability derivatives, nine response "// &
                       "solves and no finite difference")

      ! First order. `U^k` is the virtual-occupied block of the density's
      ! derivative per spin, so the total is twice it, symmetrised.
      call cphf_solve(mol, orbitals, orbital_energies, n_occ, perturbations=dip, &
                      response=u1, error=error, tol=RESPONSE_TOL)
      if (error%has_error()) return

      ! Everything below is assembled in the orbital basis and carried back.
      allocate (dm0(n_mo, n_mo), f0(n_mo, n_mo))
      dm0 = 0.0_dp
      f0 = 0.0_dp
      do i = 1, n_occ
         dm0(i, i) = 2.0_dp
      end do
      do i = 1, n_mo
         f0(i, i) = orbital_energies(i)
      end do
      allocate (d0(n_ao, n_ao))
      call to_ao(dm0, d0)

      allocate (dm1(n_mo, n_mo, 3), d1(n_ao, n_ao, 3), f1(n_mo, n_mo, 3))
      dm1 = 0.0_dp
      do k = 1, 3
         dm1(n_occ + 1:, :n_occ, k) = 2.0_dp*u1(:, :, k)
         dm1(:n_occ, n_occ + 1:, k) = 2.0_dp*transpose(u1(:, :, k))
         call to_ao(dm1(:, :, k), d1(:, :, k))
      end do
      do k = 1, 3
         do l = 1, 3
            alpha(k, l) = -sum(d1(:, :, l)*dip(:, :, k))
         end do
      end do

      ! The first-order Fock matrices: the field's own operator and the mean
      ! field of the density it induces.
      scratch = d1
      call response_mean_field(mol, scratch, zero_h, g, error, bounds=bounds)
      if (error%has_error()) return
      do k = 1, 3
         g(:, :, k) = g(:, :, k) + dip(:, :, k)
         call to_mo(g(:, :, k), f1(:, :, k))
      end do
      deallocate (g, scratch)

      ! Second order. Idempotency fixes the occupied-occupied and
      ! virtual-virtual blocks of the density from the first-order ones;
      ! only the block between them is solved for.
      allocate (dm2(n_mo, n_mo, 6), d2(n_ao, n_ao, 6), f2(n_mo, n_mo, 6))
      allocate (trial(n_mo, n_mo))
      dm2 = 0.0_dp
      do c = 1, 6
         k = PAIRS(1, c)
         l = PAIRS(2, c)
         call pic_gemm(dm1(:, :, k), dm1(:, :, l), trial)
         trial = 0.5_dp*(trial + transpose(trial))
         dm2(:n_occ, :n_occ, c) = -trial(:n_occ, :n_occ)
         dm2(n_occ + 1:, n_occ + 1:, c) = trial(n_occ + 1:, n_occ + 1:)
         call to_ao(dm2(:, :, c), d2(:, :, c))
      end do

      scratch = d2
      call response_mean_field(mol, scratch, zero_h, g, error, bounds=bounds)
      if (error%has_error()) return
      allocate (rhs(n_vir, n_occ, 6))
      do c = 1, 6
         k = PAIRS(1, c)
         l = PAIRS(2, c)
         call to_mo(g(:, :, c), f2(:, :, c))
         ! The commutators `[F^k, P^l] + [F^l, P^k]` in the virtual-occupied
         ! block, `P` being the density per spin.
         rhs(:, :, c) = f2(n_occ + 1:, :n_occ, c) &
                        + matmul(f1(n_occ + 1:, n_occ + 1:, k), u1(:, :, l)) &
                        - matmul(u1(:, :, l), f1(:n_occ, :n_occ, k)) &
                        + matmul(f1(n_occ + 1:, n_occ + 1:, l), u1(:, :, k)) &
                        - matmul(u1(:, :, k), f1(:n_occ, :n_occ, l))
      end do
      deallocate (g, scratch)
      call cphf_solve(mol, orbitals, orbital_energies, n_occ, response=u2, error=error, &
                      tol=RESPONSE_TOL, mo_rhs=rhs)
      if (error%has_error()) then
         call error%add_context("Raman: the second-order field response")
         return
      end if

      ! The solved block, and the mean field of it on top of the fixed
      ! blocks' one.
      allocate (scratch(n_ao, n_ao, 6))
      do c = 1, 6
         trial = 0.0_dp
         trial(n_occ + 1:, :n_occ) = 2.0_dp*u2(:, :, c)
         trial(:n_occ, n_occ + 1:) = 2.0_dp*transpose(u2(:, :, c))
         dm2(:, :, c) = dm2(:, :, c) + trial
         call to_ao(trial, scratch(:, :, c))
         d2(:, :, c) = d2(:, :, c) + scratch(:, :, c)
      end do
      call response_mean_field(mol, scratch, zero_h, g, error, bounds=bounds)
      if (error%has_error()) return
      do c = 1, 6
         call to_mo(g(:, :, c), trial)
         f2(:, :, c) = f2(:, :, c) + trial
      end do
      deallocate (g, scratch)

      ! The energy-weighted density is `(1/2) D F D`, differentiated by the
      ! product rule over its three factors.
      allocate (w2(n_ao, n_ao, 6), wm(n_mo, n_mo))
      do c = 1, 6
         k = PAIRS(1, c)
         l = PAIRS(2, c)
         wm = triple(dm2(:, :, c), f0, dm0) + triple(dm0, f0, dm2(:, :, c)) &
              + triple(dm0, f2(:, :, c), dm0) &
              + triple(dm1(:, :, k), f1(:, :, l), dm0) + triple(dm1(:, :, l), f1(:, :, k), dm0) &
              + triple(dm0, f1(:, :, k), dm1(:, :, l)) + triple(dm0, f1(:, :, l), dm1(:, :, k)) &
              + triple(dm1(:, :, k), f0, dm1(:, :, l)) + triple(dm1(:, :, l), f0, dm1(:, :, k))
         wm = 0.5_dp*wm
         call to_ao(wm, w2(:, :, c))
      end do
      deallocate (dm0, dm1, dm2, f0, f1, f2, wm, trial, u1, u2, rhs)

      ! The integral derivatives, in libcint's convention: the nabla is on
      ! the bra, and the derivative with respect to the bra's atom is minus
      ! that.
      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)
      allocate (offsets(natm), counts(natm))
      call atom_ao_blocks(mol, offsets, counts)

      ! The two-electron term is quadratic in the density it is contracted
      ! with, and the contraction available is the gradient's own, at one
      ! density. A product of two different densities is the difference of
      ! that at their sum and at each.
      allocate (q1(3, natm, 3), qpair(3, natm), plus(3, natm), minus(3, natm))
      q1 = 0.0_dp
      do k = 1, 3
         call two_electron_gradient(mol, d1(:, :, k), q1(:, :, k), error)
         if (error%has_error()) return
      end do

      allocate (dalpha(6, 3*natm), second(3, natm), trial(n_ao, n_ao))
      do c = 1, 6
         k = PAIRS(1, c)
         l = PAIRS(2, c)
         second = 0.0_dp

         call hellmann_feynman_gradient(mol, d2(:, :, c), second)
         do iatom = 1, natm
            if (counts(iatom) == 0) cycle
            p0 = offsets(iatom) + 1
            p1 = offsets(iatom) + counts(iatom)
            do comp = 1, 3
               second(comp, iatom) = second(comp, iatom) &
                                     + 2.0_dp*sum(h1(p0:p1, :, comp)*d2(p0:p1, :, c)) &
                                     - 2.0_dp*sum(s1(p0:p1, :, comp)*w2(p0:p1, :, c))
            end do
         end do

         ! `Tr(D^kl G^x[D])`
         plus = 0.0_dp
         minus = 0.0_dp
         trial = d0 + d2(:, :, c)
         call two_electron_gradient(mol, trial, plus, error)
         if (error%has_error()) return
         trial = d0 - d2(:, :, c)
         call two_electron_gradient(mol, trial, minus, error)
         if (error%has_error()) return
         second = second + 0.5_dp*(plus - minus)

         ! `Tr(D^k G^x[D^l])`
         if (k == l) then
            second = second + 2.0_dp*q1(:, :, k)
         else
            qpair = 0.0_dp
            trial = d1(:, :, k) + d1(:, :, l)
            call two_electron_gradient(mol, trial, qpair, error)
            if (error%has_error()) return
            second = second + qpair - q1(:, :, k) - q1(:, :, l)
         end if

         ! The field term's own derivative, at the first-order densities.
         do iatom = 1, natm
            do comp = 1, 3
               second(comp, iatom) = second(comp, iatom) &
                                     + sum(d1(:, :, k)*ddip(:, :, l, comp, iatom)) &
                                     + sum(d1(:, :, l)*ddip(:, :, k, comp, iatom))
               dalpha(c, 3*(iatom - 1) + comp) = -second(comp, iatom)
            end do
         end do
      end do

   contains

      subroutine to_ao(mo, ao)
         !! `C mo C^T`
         real(dp), intent(in) :: mo(:, :)
         real(dp), intent(out) :: ao(:, :)
         real(dp), allocatable :: half(:, :)

         allocate (half(n_ao, n_mo))
         call pic_gemm(orbitals, mo, half)
         call pic_gemm(half, orbitals, ao, transb="T")
      end subroutine to_ao

      subroutine to_mo(ao, mo)
         !! `C^T ao C`
         real(dp), intent(in) :: ao(:, :)
         real(dp), intent(out) :: mo(:, :)
         real(dp), allocatable :: half(:, :)

         allocate (half(n_ao, n_mo))
         call pic_gemm(ao, orbitals, half)
         call pic_gemm(orbitals, half, mo, transa="T")
      end subroutine to_mo

      function triple(x, y, z) result(xyz)
         !! The product of three square matrices
         real(dp), intent(in) :: x(:, :), y(:, :), z(:, :)
         real(dp) :: xyz(size(x, 1), size(x, 1))
         real(dp), allocatable :: xy(:, :)

         allocate (xy(size(x, 1), size(x, 1)))
         call pic_gemm(x, y, xy)
         call pic_gemm(xy, z, xyz)
      end function triple
   end subroutine polarizability_derivatives_analytic

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
         call run_czt_rhf(mol, nelec, settings%max_iter, SCF_ENERGY_TOL, SCF_DENSITY_TOL, &
                          .false., scf, error, guess=SCF_GUESS_SAD, &
                          guess_density=guess_density, xc=xc, &
                          level_shift=settings%level_shift, &
                          linear_dependence=settings%linear_dependence, &
                          accelerator=accelerator, diis_vectors=diis_vectors, &
                          incremental_fock=settings%incremental_fock, &
                          grad_tol=convergence%gradient_tolerance, &
                          convergence=convergence)
      else
         call run_czt_rhf(mol, nelec, settings%max_iter, SCF_ENERGY_TOL, SCF_DENSITY_TOL, &
                          .false., scf, error, guess=SCF_GUESS_SAD, &
                          guess_density=guess_density, &
                          level_shift=settings%level_shift, &
                          linear_dependence=settings%linear_dependence, &
                          accelerator=accelerator, diis_vectors=diis_vectors, &
                          incremental_fock=settings%incremental_fock, &
                          grad_tol=convergence%gradient_tolerance, &
                          convergence=convergence)
      end if
      if (.not. error%has_error() .and. .not. scf%converged) then
         call error%set(ERROR_VALIDATION, "the reference did not converge at a "// &
                        "displaced geometry")
      end if
      if (.not. error%has_error()) then
         if (kohn_sham) then
            call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                       scf%n_occupied, alpha, error, tol=RESPONSE_TOL, &
                                       xc=xc, density=scf%density)
         else
            call static_polarizability(mol, scf%orbitals, scf%orbital_energies, &
                                       scf%n_occupied, alpha, error, tol=RESPONSE_TOL)
         end if
      end if
      if (kohn_sham) call xc%destroy()
      call mol%destroy()
   end subroutine polarizability_at

end module mqc_czt_raman
