!! Basis function values on a grid, for DFT
module mqc_libcint_ao
   !! Evaluates every basis function of a molecule at a set of points.
   !!
   !! **libcint does not do this.** It computes integrals; there is no `GTOval` in
   !! it, and the function of that name belongs to PySCF's own `libcgto`. So this
   !! is the one piece of the DFT path that had to be written rather than called,
   !! and the whole difficulty is matching libcint's conventions exactly -- a basis
   !! function evaluated under a different-but-valid convention is a different
   !! basis, and the resulting SCF converges tidily onto a wrong number.
   !!
   !! Three conventions have to line up, and each is a way to be wrong quietly:
   !!
   !!   * **Normalisation is already in `env`.** `molecule_build` folds
   !!     `libcint_gto_norm` into every contraction coefficient, because libcint
   !!     expects that and does not apply it itself. Nothing here may normalise
   !!     again; doing so gives an overlap diagonal that is not one.
   !!   * **Cartesian component order** is libcint's: for l, the components run
   !!     `x^(l-i) y^(i-j) z^j` over `i = 0..l`, `j = 0..i`, which for d is
   !!     xx, xy, xz, yy, yz, zz.
   !!   * **The spherical transform must be libcint's**, coefficients and all. They
   !!     are in `mqc_libcint_ao_data`, transcribed from its own table, with the
   !!     s and p normalisation it keeps outside that table applied here.
   !!
   !! Which convention a molecule is in comes from `mol%cartesian`, the same flag
   !! `shell_dim` routes on, so the AO values cannot disagree with the AO count.
   !!
   !! Points are taken in blocks. A medium grid on a modest molecule is tens of
   !! thousands of points, and `n_points` by `n_ao` held whole is hundreds of
   !! megabytes at a real basis size -- the shape of mistake the coupled-cluster
   !! work already paid for once.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   use mqc_libcint_ao_data, only: C2S_LMAX, c2s_block, common_fac_sp
   use libcint_fortran, only: LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF
   implicit none
   private

   public :: eval_ao_block
   public :: eval_rho
   public :: max_ao_l

   !> Points per pass when a caller asks for a whole grid at once.
   integer, parameter, public :: AO_POINT_BLOCK = 4096

contains

   pure function max_ao_l(mol) result(l_max)
      !! Highest angular momentum in this molecule's basis
      type(libcint_molecule_t), intent(in) :: mol
      integer :: l_max

      integer :: ish

      l_max = 0
      do ish = 1, mol%nbas
         l_max = max(l_max, mol%bas(LIBCINT_ANG_OF, ish))
      end do
   end function max_ao_l

   subroutine eval_ao_block(mol, coords, ao, error, grad)
      !! chi_mu(r) for every basis function at every supplied point
      !!
      !! `ao` comes back as (n_points, n_ao), which is the orientation the density
      !! assembly wants: rho at a point is a row contracted against the density
      !! matrix, so points vary slowest and a gemm over them reads contiguously.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)        !! (3, n_points), Bohr
      real(dp), allocatable, intent(out) :: ao(:, :)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: grad(:, :, :)
         !! d chi / d r, as (n_points, n_ao, 3). Asked for by GGA and above.
         !!
         !! Produced here rather than by a second routine because it shares
         !! everything expensive: the same exponentials, the same angular
         !! components, the same spherical transform. A separate `eval_ao_deriv1`
         !! would duplicate all three and have to be kept in step with them.

      integer :: n_points, ish, l, nprim, nctr, ic, ip, ig, off_ao
      integer :: n_cart, n_sph, n_here, comp, i, j, iatom
      integer :: exp_ptr, coef_ptr
      real(dp) :: dx, dy, dz, r2, radial, fac, dradial
      real(dp), allocatable :: cart(:), sph(:), trans(:, :)
      real(dp), allocatable :: dcart(:, :), dsph(:, :)
      logical :: want_grad
      integer :: id

      n_points = size(coords, 2)

      if (max_ao_l(mol) > C2S_LMAX) then
         call error%set(ERROR_VALIDATION, "AO evaluation: this basis has an angular "// &
                        "momentum above the transform table's range. Extend "// &
                        "mqc_libcint_ao_data rather than letting it through -- the "// &
                        "functions would otherwise be silently wrong.")
         return
      end if

      want_grad = present(grad)
      allocate (ao(n_points, mol%nao))
      ao = 0.0_dp
      if (want_grad) then
         allocate (grad(n_points, mol%nao, 3))
         grad = 0.0_dp
      end if

      !$omp parallel default(none) shared(mol, coords, ao, grad, n_points, want_grad) &
      !$omp    private(ish, l, nprim, nctr, ic, ip, ig, off_ao, n_cart, n_sph, &
      !$omp            n_here, comp, i, j, iatom, exp_ptr, coef_ptr, id, &
      !$omp            dx, dy, dz, r2, radial, dradial, fac, cart, sph, trans, &
      !$omp            dcart, dsph)
      allocate (cart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2), sph(2*C2S_LMAX + 1))
      allocate (trans(2*C2S_LMAX + 1, (C2S_LMAX + 1)*(C2S_LMAX + 2)/2))
      allocate (dcart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2, 3), dsph(2*C2S_LMAX + 1, 3))

      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         nprim = mol%bas(LIBCINT_NPRIM_OF, ish)
         nctr = mol%bas(LIBCINT_NCTR_OF, ish)
         exp_ptr = mol%bas(LIBCINT_PTR_EXP, ish)
         coef_ptr = mol%bas(LIBCINT_PTR_COEFF, ish)
         iatom = mol%bas(LIBCINT_ATOM_OF, ish) + 1

         n_cart = (l + 1)*(l + 2)/2
         n_sph = 2*l + 1
         n_here = shell_dim(mol%cartesian, ish - 1, mol%bas)/nctr

         ! The factor libcint keeps outside its transform table, and it applies to
         ! *both* angular conventions -- `cint1e.c` multiplies every integral by
         ! `CINTcommon_fac_sp` of each shell's l without consulting cart-vs-sph.
         ! Gating it on the spherical path is wrong and shows up only in a
         ! Cartesian basis, where s and p come out too large by 1/0.2820948 and
         ! 1/0.4886025.
         fac = common_fac_sp(l)
         if (.not. mol%cartesian) call c2s_block(l, trans)

         !$omp do schedule(static)
         do ig = 1, n_points
            dx = coords(1, ig) - mol%coords(1, iatom)
            dy = coords(2, ig) - mol%coords(2, iatom)
            dz = coords(3, ig) - mol%coords(3, iatom)
            r2 = dx*dx + dy*dy + dz*dz

            ! Angular part, in libcint's Cartesian order, and its gradient. The
            ! power rule term vanishes when the exponent is zero, which has to be
            ! a branch rather than an evaluation: pow(x, -1) is not what it means.
            comp = 0
            do i = 0, l
               do j = 0, i
                  comp = comp + 1
                  cart(comp) = pow(dx, l - i)*pow(dy, i - j)*pow(dz, j)
                  if (want_grad) then
                     dcart(comp, :) = 0.0_dp
                     if (l - i > 0) dcart(comp, 1) = real(l - i, dp) &
                                                     *pow(dx, l - i - 1)*pow(dy, i - j)*pow(dz, j)
                     if (i - j > 0) dcart(comp, 2) = real(i - j, dp) &
                                                     *pow(dx, l - i)*pow(dy, i - j - 1)*pow(dz, j)
                     if (j > 0) dcart(comp, 3) = real(j, dp) &
                                                 *pow(dx, l - i)*pow(dy, i - j)*pow(dz, j - 1)
                  end if
               end do
            end do

            ! The transform is linear and does not depend on position, so it
            ! applies to the gradient exactly as it does to the value.
            if (.not. mol%cartesian) then
               do comp = 1, n_sph
                  sph(comp) = fac*sum(trans(comp, 1:n_cart)*cart(1:n_cart))
                  if (want_grad) then
                     do id = 1, 3
                        dsph(comp, id) = fac*sum(trans(comp, 1:n_cart)*dcart(1:n_cart, id))
                     end do
                  end if
               end do
            end if

            ! One contraction column at a time; libcint lays a shell's functions
            ! out with the contraction index outermost, and `molecule_build`
            ! preserved that, so the AO offset advances by n_here per column.
            do ic = 1, nctr
               radial = 0.0_dp
               dradial = 0.0_dp
               do ip = 1, nprim
                  ! d/dr_k exp(-a r^2) = -2 a r_k exp(-a r^2), so the -2a is
                  ! accumulated here and the r_k factor applied per component.
                  radial = radial + mol%env(coef_ptr + (ic - 1)*nprim + ip) &
                           *exp(-mol%env(exp_ptr + ip)*r2)
                  if (want_grad) then
                     dradial = dradial - 2.0_dp*mol%env(exp_ptr + ip) &
                               *mol%env(coef_ptr + (ic - 1)*nprim + ip) &
                               *exp(-mol%env(exp_ptr + ip)*r2)
                  end if
               end do
               off_ao = mol%shell_offset(ish) + (ic - 1)*n_here
               if (mol%cartesian) then
                  do comp = 1, n_cart
                     ao(ig, off_ao + comp) = fac*radial*cart(comp)
                  end do
                  if (want_grad) then
                     do comp = 1, n_cart
                        grad(ig, off_ao + comp, 1) = fac*(dcart(comp, 1)*radial &
                                                          + cart(comp)*dradial*dx)
                        grad(ig, off_ao + comp, 2) = fac*(dcart(comp, 2)*radial &
                                                          + cart(comp)*dradial*dy)
                        grad(ig, off_ao + comp, 3) = fac*(dcart(comp, 3)*radial &
                                                          + cart(comp)*dradial*dz)
                     end do
                  end if
               else
                  do comp = 1, n_sph
                     ao(ig, off_ao + comp) = radial*sph(comp)
                  end do
                  if (want_grad) then
                     do comp = 1, n_sph
                        grad(ig, off_ao + comp, 1) = dsph(comp, 1)*radial &
                                                     + sph(comp)*dradial*dx
                        grad(ig, off_ao + comp, 2) = dsph(comp, 2)*radial &
                                                     + sph(comp)*dradial*dy
                        grad(ig, off_ao + comp, 3) = dsph(comp, 3)*radial &
                                                     + sph(comp)*dradial*dz
                     end do
                  end if
               end if
            end do
         end do
         !$omp end do
      end do

      deallocate (cart, sph, trans, dcart, dsph)
      !$omp end parallel
   end subroutine eval_ao_block

   subroutine eval_rho(ao, density, rho, ao_grad, rho_grad, tau)
      !! The electron density at every point, from the density matrix
      !!
      !!     rho(r) = sum_uv D_uv chi_u(r) chi_v(r)
      !!
      !! Through a gemm rather than a double loop over basis functions: form
      !! X = chi D once, then rho at a point is the row-wise dot of X with chi.
      !! The alternative is n_points n_ao^2 scalar work with the density matrix
      !! read in the worst possible order, which is the mistake the coupled-cluster
      !! ladder already paid for.
      !!
      !! `density` is whatever convention the caller's SCF uses -- the restricted
      !! path's D already carries two electrons per orbital, so this returns the
      !! total density for it and a spin density for an unrestricted D. Nothing
      !! here applies a factor, because doing so would be right for exactly one of
      !! those and silently wrong for the other.
      real(dp), intent(in) :: ao(:, :)        !! (n_points, n_ao)
      real(dp), intent(in) :: density(:, :)   !! (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: rho(:)
      real(dp), intent(in), optional :: ao_grad(:, :, :)
      real(dp), allocatable, intent(out), optional :: rho_grad(:, :)
         !! d rho / d r, as (n_points, 3). Both optionals or neither.
         !!
         !!     grad rho = 2 sum_uv D_uv chi_v grad chi_u
         !!
         !! The factor of two is the symmetry of D, not a spin factor, and it is
         !! the term that a GGA gets wrong by a few millihartree while converging
         !! perfectly.

      real(dp), allocatable, intent(out), optional :: tau(:)
         !! The kinetic energy density, for meta-GGA. Needs `ao_grad`.
         !!
         !!     tau(r) = 1/2 sum_d sum_uv D_uv d_d chi_u d_d chi_v
         !!
         !! The half is libxc's convention and PySCF's, checked against
         !! `eval_rho(..., xctype='MGGA')` rather than assumed -- the same quantity
         !! is defined without it elsewhere in the literature, and a factor of two
         !! in tau converges to a wrong energy like everything else here.

      real(dp), allocatable :: work(:, :), gwork(:, :)
      integer :: n_points, n_ao, ig, mu, id

      n_points = size(ao, 1)
      n_ao = size(ao, 2)

      allocate (work(n_points, n_ao), rho(n_points))
      call pic_gemm(ao, density, work)

      rho = 0.0_dp
      !$omp parallel do default(none) shared(rho, work, ao, n_points, n_ao) &
      !$omp    private(ig, mu) schedule(static)
      do ig = 1, n_points
         do mu = 1, n_ao
            rho(ig) = rho(ig) + work(ig, mu)*ao(ig, mu)
         end do
      end do
      !$omp end parallel do

      if (present(tau)) then
         allocate (tau(n_points))
         tau = 0.0_dp
         if (.not. present(ao_grad)) then
            tau = huge(1.0_dp)
         else
            allocate (gwork(n_points, n_ao))
            do id = 1, 3
               call pic_gemm(ao_grad(:, :, id), density, gwork)
               !$omp parallel do default(none) &
               !$omp    shared(tau, gwork, ao_grad, n_points, n_ao, id) &
               !$omp    private(ig, mu) schedule(static)
               do ig = 1, n_points
                  do mu = 1, n_ao
                     tau(ig) = tau(ig) + 0.5_dp*gwork(ig, mu)*ao_grad(ig, mu, id)
                  end do
               end do
               !$omp end parallel do
            end do
            deallocate (gwork)
         end if
      end if

      if (present(rho_grad)) then
         allocate (rho_grad(n_points, 3))
         rho_grad = 0.0_dp
         if (.not. present(ao_grad)) then
            ! Nothing to build it from; leaving zeros would be a silently
            ! LDA-shaped answer under a GGA's name.
            rho_grad = huge(1.0_dp)
         else
            !$omp parallel do default(none) &
            !$omp    shared(rho_grad, work, ao_grad, n_points, n_ao) &
            !$omp    private(ig, mu, id) schedule(static)
            do ig = 1, n_points
               do id = 1, 3
                  do mu = 1, n_ao
                     rho_grad(ig, id) = rho_grad(ig, id) &
                                        + 2.0_dp*work(ig, mu)*ao_grad(ig, mu, id)
                  end do
               end do
            end do
            !$omp end parallel do
         end if
      end if

      deallocate (work)
   end subroutine eval_rho

   pure function pow(x, n) result(v)
      !! x**n for small non-negative n, without the intrinsic's generality
      !!
      !! `x**0` must be one even when x is zero, which is exactly what happens at
      !! a grid point sitting on a nucleus -- and there are such points.
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: v

      integer :: k

      v = 1.0_dp
      do k = 1, n
         v = v*x
      end do
   end function pow

end module mqc_libcint_ao
