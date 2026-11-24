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
   use pic_io, only: to_char
   use pic_timer, only: timer_type
   use pic_logger, only: logger => global_logger
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, PROPERTY_RESPONSE_TOL, &
                                       DIFFERENCED_SCF_ENERGY_TOL, &
                                       DIFFERENCED_SCF_DENSITY_TOL
   use mqc_symmetric_tensor, only: pack_symmetric, N_SYMMETRIC, SYMMETRIC_ROW, &
                                   SYMMETRIC_COL
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf, SCF_GUESS_SAD
   use pic_blas_interfaces, only: pic_gemm
   use mqc_czt_integrals, only: atom_ao_blocks
   use mqc_czt_cphf, only: static_polarizability, cphf_solve, build_hessian_mo, &
                           mo_transform_fits, factorize_response
   use mqc_memory, only: memory_budget
   use omp_lib, only: omp_get_max_threads
   use pic_lapack_interfaces, only: pic_getrs
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_response_product, only: response_mean_field
   use mqc_czt_gradient, only: czt_scf_gradient, one_electron_deriv, &
                               two_electron_deriv_many, hellmann_feynman_gradient, &
                               xc_potential_gradient, xc_kernel_gradient, &
                               DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_czt_multipole, only: multipole_matrices, dipole_integral_derivatives
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_kernel2_apply
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_scf_convergence, only: scf_convergence_t
   implicit none
   private

   public :: polarizability_derivatives
   public :: polarizability_derivatives_field
   public :: polarizability_derivatives_analytic

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
   integer, parameter :: N_STACK = 1 + 3 + N_SYMMETRIC
      !! Densities whose derivative potentials the analytic route needs: the
      !! reference's, the three first-order ones and the six second-order ones
   real(dp), parameter :: STACK_BUDGET_SHARE = 0.5_dp
      !! Share of the memory budget the derivative pass over the stack may
      !! take; what it does not fit in goes through in pieces.
   real(dp), parameter :: STACK_BLIND_BUDGET = 4.0e9_dp
      !! Bytes that pass plans on where the machine reports nothing

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
                                               alpha, dalpha, error, xc, field_strength)
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
      real(dp), intent(in), optional :: field_strength
         !! The field the differences are taken at, atomic units.
         !! `RAMAN_FIELD` when absent. The error falls as its square, so two
         !! strengths extrapolate.

      ! The off-diagonal components, in the packed order xy, xz, yz.
      integer, parameter :: OFF_DIAGONAL(2, 3) = reshape([1, 2, 1, 3, 2, 3], [2, 3])
      real(dp), allocatable :: dip(:, :, :), ddip(:, :, :, :, :)
      real(dp), allocatable :: g0(:, :), g(:, :, :), second(:, :, :)
      real(dp) :: origin(3), fields(3, N_FIELDS)
      type(scf_convergence_t) :: tight
      type(rhf_result_t) :: scf
      integer :: natm, k, l, ipair, ifield, iatom, comp, col
      logical :: kohn_sham
      real(dp) :: step

      alpha = 0.0_dp
      if (error%has_error()) return
      step = RAMAN_FIELD
      if (present(field_strength)) step = field_strength
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
         fields(k, 2*k - 1) = step
         fields(k, 2*k) = -step
      end do
      do ipair = 1, 3
         fields(OFF_DIAGONAL(:, ipair), 5 + 2*ipair) = step
         fields(OFF_DIAGONAL(:, ipair), 6 + 2*ipair) = -step
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
         second(:, :, k) = (g(:, :, 2*k - 1) + g(:, :, 2*k) - 2.0_dp*g0)/step**2
      end do
      do ipair = 1, 3
         k = OFF_DIAGONAL(1, ipair)
         l = OFF_DIAGONAL(2, ipair)
         second(:, :, 3 + ipair) = 0.5_dp*((g(:, :, 5 + 2*ipair) + g(:, :, 6 + 2*ipair) &
                                            - 2.0_dp*g0)/step**2 &
                                           - second(:, :, k) - second(:, :, l))
      end do

      ! Components xx, yy, zz, xy, xz, yz, which is the order of `OFF_DIAGONAL`.
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

   subroutine polarizability_derivatives_analytic(mol, orbitals, orbital_energies, n_occ, &
                                                  alpha, dalpha, error, integral_direct, &
                                                  stored_response, xc)
      !! `d alpha / d x` with no finite difference, for a closed shell under
      !! Hartree-Fock or a Kohn-Sham LDA, GGA or global hybrid
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
      !! and one pass over the differentiated two-electron integrals.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)          !! Canonical, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), intent(out) :: alpha(3, 3)   !! Bohr^3
      real(dp), allocatable, intent(out) :: dalpha(:, :)
         !! (6, 3*n_atoms) in Bohr^2, laid out as `polarizability_derivatives`
         !! lays it out
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: integral_direct
         !! Solve the nine response equations by iteration, recomputing the
         !! integrals on every pass, even where the response matrix would fit
         !! in memory. Off when absent; the two give the same derivatives.
      logical, intent(out), optional :: stored_response
         !! Whether the solves came from the factorized response matrix
      type(xc_context_t), intent(inout), optional, target :: xc
         !! The reference's functional. Absent or inactive is Hartree-Fock.
         !! A meta-GGA, a range-separated hybrid and non-local correlation
         !! are refused: their third derivatives are not built.

      real(dp), allocatable :: dip(:, :, :), ddip(:, :, :, :, :)
      real(dp), allocatable :: aplus(:, :), aminus(:, :), no_eri(:, :, :, :), gaps(:, :)
      integer, allocatable :: ipiv(:)
      integer :: info
      logical :: stored
      real(dp), allocatable :: u1(:, :, :), u2(:, :, :), rhs(:, :, :)
      real(dp), allocatable :: d0(:, :), d1(:, :, :), d2(:, :, :), w2(:, :, :)
      real(dp), allocatable :: dm1(:, :, :), dm2(:, :, :), f1(:, :, :), f2(:, :, :)
      real(dp), allocatable :: dm0(:, :), f0(:, :), wm(:, :)
      real(dp), allocatable :: g(:, :, :), scratch(:, :, :), bounds(:, :), zero_h(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :)
      real(dp), allocatable :: stack(:, :, :), vhf(:, :, :, :), vhf_part(:, :, :, :)
      real(dp) :: per_set
      integer :: n_chunk, first_set, last_set
      real(dp) :: exchange
      logical :: kohn_sham
      type(xc_context_t), pointer :: functional
      real(dp), allocatable :: second(:, :), trial(:, :)
      integer, allocatable :: offsets(:), counts(:)
      real(dp) :: origin(3)
      integer :: n_ao, n_mo, n_vir, natm, k, l, c, i, iatom, comp, p0, p1
      type(timer_type) :: clock

      alpha = 0.0_dp
      if (present(stored_response)) stored_response = .false.
      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_vir = n_mo - n_occ
      natm = mol%natm

      ! Under a functional the mean field keeps a fraction of exchange and
      ! gains the kernel; the energy gains terms in the functional's second
      ! and third derivatives.
      kohn_sham = .false.
      if (present(xc)) kohn_sham = xc%active
      ! What goes down to the response routines: the context when it carries
      ! a functional, and nothing otherwise. An inactive context is present,
      ! and one that was never created says its exchange fraction is zero.
      functional => null()
      exchange = 1.0_dp
      if (kohn_sham) then
         functional => xc
         if (xc%any_mgga .or. xc%range_separated .or. xc%nlc_b /= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "Raman: the analytic polarizability "// &
                           "derivatives are built for LDA, GGA and global hybrid "// &
                           "functionals; a meta-GGA, a range-separated hybrid and "// &
                           "non-local correlation need third derivatives that are not")
            return
         end if
         exchange = xc%exx_fraction
      end if

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
      call clock%start()

      ! Nine right-hand sides of one operator, `(A+B)`. Where that matrix
      ! fits it is built once by transforming the integrals and factorized,
      ! and each solve is a back-substitution; otherwise each is iterated,
      ! with a pass over the integrals per iteration.
      stored = mo_transform_fits(n_ao, n_occ, n_vir, .true.)
      if (present(integral_direct)) stored = stored .and. .not. integral_direct
      ! The transformed matrix is the Hartree-Fock one: it has no kernel.
      if (kohn_sham) stored = .false.
      if (present(stored_response)) stored_response = stored
      if (stored) then
         allocate (gaps(n_vir, n_occ), no_eri(0, 0, 0, 0))
         do i = 1, n_occ
            gaps(:, i) = orbital_energies(n_occ + 1:) - orbital_energies(i)
         end do
         ! What the iterated solve refuses, this one has to refuse too: a
         ! factorization goes through on an indefinite matrix without a word.
         if (minval(gaps) <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "Raman: the reference has no gap between "// &
                           "its occupied and virtual orbitals")
            return
         end if
         call build_hessian_mo(mol, no_eri, orbitals(:, :n_occ), orbitals(:, n_occ + 1:), &
                               gaps, aplus, aminus, error)
         if (error%has_error()) return
         deallocate (aminus, gaps, no_eri)
         allocate (ipiv(n_vir*n_occ))
         call factorize_response(n_vir*n_occ, aplus, ipiv, info)
         if (info /= 0) then
            call error%set(ERROR_VALIDATION, "Raman: the response matrix is singular")
            return
         end if
      end if

      allocate (d0(n_ao, n_ao))
      d0 = 2.0_dp*matmul(orbitals(:, :n_occ), transpose(orbitals(:, :n_occ)))

      ! First order. `U^k` is the virtual-occupied block of the density's
      ! derivative per spin, so the total is twice it, symmetrised.
      if (stored) then
         allocate (rhs(n_vir, n_occ, 3), trial(n_mo, n_mo))
         do k = 1, 3
            call to_mo(dip(:, :, k), trial)
            rhs(:, :, k) = trial(n_occ + 1:, :n_occ)
         end do
         call solve_stored(rhs, u1)
         deallocate (rhs, trial)
      else
         call cphf_solve(mol, orbitals, orbital_energies, n_occ, perturbations=dip, &
                         response=u1, error=error, tol=PROPERTY_RESPONSE_TOL, &
                         xc=functional, density=d0)
      end if
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
      call response_mean_field(mol, scratch, zero_h, g, error, bounds=bounds, &
                               xc=functional, reference=d0)
      if (error%has_error()) return
      do k = 1, 3
         g(:, :, k) = g(:, :, k) + dip(:, :, k)
         call to_mo(g(:, :, k), f1(:, :, k))
      end do
      deallocate (g, scratch)

      ! Second order. Idempotency fixes the occupied-occupied and
      ! virtual-virtual blocks of the density from the first-order ones;
      ! only the block between them is solved for.
      allocate (dm2(n_mo, n_mo, N_SYMMETRIC), d2(n_ao, n_ao, N_SYMMETRIC), f2(n_mo, n_mo, N_SYMMETRIC))
      allocate (trial(n_mo, n_mo))
      dm2 = 0.0_dp
      do c = 1, N_SYMMETRIC
         k = SYMMETRIC_ROW(c)
         l = SYMMETRIC_COL(c)
         call pic_gemm(dm1(:, :, k), dm1(:, :, l), trial)
         trial = 0.5_dp*(trial + transpose(trial))
         dm2(:n_occ, :n_occ, c) = -trial(:n_occ, :n_occ)
         dm2(n_occ + 1:, n_occ + 1:, c) = trial(n_occ + 1:, n_occ + 1:)
         call to_ao(dm2(:, :, c), d2(:, :, c))
      end do

      scratch = d2
      call response_mean_field(mol, scratch, zero_h, g, error, bounds=bounds, &
                               xc=functional, reference=d0)
      if (error%has_error()) return
      allocate (rhs(n_vir, n_occ, N_SYMMETRIC))
      do c = 1, N_SYMMETRIC
         k = SYMMETRIC_ROW(c)
         l = SYMMETRIC_COL(c)
         ! The potential is not linear in the density: its second variation
         ! along the two first-order densities is part of `F^kl` whatever
         ! `D^kl` turns out to be.
         if (kohn_sham) then
            call xc_kernel2_apply(functional, mol, d0, d1(:, :, k), d1(:, :, l), g(:, :, c), error)
            if (error%has_error()) return
         end if
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
      if (stored) then
         call solve_stored(rhs, u2)
         ! The last use of the factors, and the largest thing held.
         if (allocated(aplus)) deallocate (aplus)
         if (allocated(ipiv)) deallocate (ipiv)
      else
         call cphf_solve(mol, orbitals, orbital_energies, n_occ, response=u2, error=error, &
                         tol=PROPERTY_RESPONSE_TOL, mo_rhs=rhs, xc=functional, density=d0)
      end if
      if (error%has_error()) then
         call error%add_context("Raman: the second-order field response")
         return
      end if

      ! The solved block, and the mean field of it on top of the fixed
      ! blocks' one.
      allocate (scratch(n_ao, n_ao, N_SYMMETRIC))
      do c = 1, N_SYMMETRIC
         trial = 0.0_dp
         trial(n_occ + 1:, :n_occ) = 2.0_dp*u2(:, :, c)
         trial(:n_occ, n_occ + 1:) = 2.0_dp*transpose(u2(:, :, c))
         dm2(:, :, c) = dm2(:, :, c) + trial
         call to_ao(trial, scratch(:, :, c))
         d2(:, :, c) = d2(:, :, c) + scratch(:, :, c)
      end do
      call response_mean_field(mol, scratch, zero_h, g, error, bounds=bounds, &
                               xc=functional, reference=d0)
      if (error%has_error()) return
      do c = 1, N_SYMMETRIC
         call to_mo(g(:, :, c), trial)
         f2(:, :, c) = f2(:, :, c) + trial
      end do
      deallocate (g, scratch)

      ! The energy-weighted density is `(1/2) D F D`, differentiated by the
      ! product rule over its three factors.
      allocate (w2(n_ao, n_ao, N_SYMMETRIC), wm(n_mo, n_mo))
      do c = 1, N_SYMMETRIC
         k = SYMMETRIC_ROW(c)
         l = SYMMETRIC_COL(c)
         wm = triple(dm2(:, :, c), f0, dm0) + triple(dm0, f0, dm2(:, :, c)) &
              + triple(dm0, f2(:, :, c), dm0) &
              + triple(dm1(:, :, k), f1(:, :, l), dm0) + triple(dm1(:, :, l), f1(:, :, k), dm0) &
              + triple(dm0, f1(:, :, k), dm1(:, :, l)) + triple(dm0, f1(:, :, l), dm1(:, :, k)) &
              + triple(dm1(:, :, k), f0, dm1(:, :, l)) + triple(dm1(:, :, l), f0, dm1(:, :, k))
         wm = 0.5_dp*wm
         call to_ao(wm, w2(:, :, c))
      end do
      deallocate (dm0, dm1, dm2, f0, f1, f2, wm, trial, u1, u2, rhs)

      call logger%info("  Raman: response solves and densities in "// &
                       to_char(clock%get_elapsed_time())//" s")
      call clock%start()

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

      ! The two-electron terms pair two different densities, `D^kl` with `D`
      ! and `D^k` with `D^l`. One pass over the differentiated integrals gives
      ! every density's derivative potential, `G^x[D]` with the nabla on the
      ! bra, and each pairing is then two traces over an atom's functions.
      allocate (stack(n_ao, n_ao, N_STACK))
      stack(:, :, 1) = d0
      stack(:, :, 2:4) = d1
      stack(:, :, 5:N_STACK) = d2
      ! The pass holds two potentials per density on every thread, so the
      ! stack goes through in as many pieces as that takes to fit.
      per_set = real(2*omp_get_max_threads() + 3, dp)*24.0_dp*real(n_ao, dp)**2
      n_chunk = int(min(real(N_STACK, dp), &
                        max(1.0_dp, memory_budget(STACK_BLIND_BUDGET, STACK_BUDGET_SHARE)/per_set)))
      allocate (vhf(n_ao, n_ao, 3, N_STACK))
      do first_set = 1, N_STACK, n_chunk
         last_set = min(first_set + n_chunk - 1, N_STACK)
         call two_electron_deriv_many(mol, stack(:, :, first_set:last_set), vhf_part, error, &
                                      exx_fraction=exchange)
         if (error%has_error()) return
         vhf(:, :, :, first_set:last_set) = vhf_part
         deallocate (vhf_part)
      end do

      allocate (dalpha(N_SYMMETRIC, 3*natm), second(3, natm))
      do c = 1, N_SYMMETRIC
         k = SYMMETRIC_ROW(c)
         l = SYMMETRIC_COL(c)
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

         ! `Tr(D^kl G^x[D])` and `Tr(D^k G^x[D^l])`
         call add_pairing(1, 4 + c)
         call add_pairing(1 + k, 1 + l)

         ! The functional's share: the potential traced with `D^kl`, and
         ! the kernel between `D^k` and `D^l`, each differentiated with the
         ! densities held and the grid moving.
         if (kohn_sham) then
            call xc_potential_gradient(functional, mol, d0, d2(:, :, c), second, error)
            call xc_kernel_gradient(functional, mol, d0, d1(:, :, k), d1(:, :, l), second, error)
            if (error%has_error()) return
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
      call logger%info("  Raman: integral derivatives contracted in "// &
                       to_char(clock%get_elapsed_time())//" s")

   contains

      subroutine solve_stored(h, u)
         !! `(A+B) U = -h` for a stack of right-hand sides, from the factors
         real(dp), intent(in) :: h(:, :, :)   !! (n_vir, n_occ, n)
         real(dp), allocatable, intent(out) :: u(:, :, :)

         real(dp), allocatable :: b(:, :), minus_h(:, :)
         integer :: status, set

         allocate (b(n_vir*n_occ, size(h, 3)), minus_h(n_vir*n_occ, size(h, 3)))
         minus_h = -reshape(h, shape(minus_h))
         b = minus_h
         call pic_getrs(aplus, ipiv, b, info=status)
         if (status /= 0) then
            call error%set(ERROR_VALIDATION, "Raman: the response equations would "// &
                           "not solve")
            return
         end if
         ! `b . (A+B)^-1 b` is positive for every `b` when the reference is a
         ! minimum. One that is not proves a saddle; all of them positive
         ! proves nothing, which is as far as the iterated solve sees too.
         do set = 1, size(h, 3)
            if (.not. any(minus_h(:, set) /= 0.0_dp)) cycle
            if (.not. sum(b(:, set)*minus_h(:, set)) > 0.0_dp) then
               call error%set(ERROR_VALIDATION, "Raman: the response matrix is not "// &
                              "positive definite; the reference SCF is a saddle "// &
                              "point, not a minimum")
               return
            end if
         end do
         allocate (u(n_vir, n_occ, size(h, 3)))
         u = reshape(b, shape(u))
      end subroutine solve_stored

      subroutine add_pairing(first, other)
         !! Add the two-electron term that pairs two densities of the stack
         !!
         !! The gradient's two-electron term at a density `D` is
         !! `2 sum_{i on A} sum_j G^x[D]_ij D_ij`. This is the part of that at
         !! `D + D'` which is linear in each: both orders, since the sum over
         !! an atom's functions is on the first index alone.
         integer, intent(in) :: first, other   !! Positions in `stack`

         integer :: atom, lo, hi, x

         do atom = 1, natm
            if (counts(atom) == 0) cycle
            lo = offsets(atom) + 1
            hi = offsets(atom) + counts(atom)
            do x = 1, 3
               second(x, atom) = second(x, atom) &
                                 + 2.0_dp*sum(vhf(lo:hi, :, x, first)*stack(lo:hi, :, other)) &
                                 + 2.0_dp*sum(vhf(lo:hi, :, x, other)*stack(lo:hi, :, first))
            end do
         end do
      end subroutine add_pairing

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
