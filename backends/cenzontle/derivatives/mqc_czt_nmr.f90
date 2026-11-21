!! Nuclear magnetic shielding tensors for a closed-shell reference
module mqc_czt_nmr
   !! The shielding of nucleus `N` is the mixed second derivative of the
   !! energy with respect to the external field `B` and that nucleus's
   !! magnetic moment `mu`,
   !!
   !!     sigma_xy = d2E / dB_x dmu_y
   !!
   !! computed over gauge-including atomic orbitals, so the result does not
   !! depend on where the origin of the vector potential is put.
   !!
   !! It comes in two parts. The diamagnetic one is an expectation value over
   !! the unperturbed density. The paramagnetic one contracts the orbitals'
   !! first-order response to the field with the nuclear-moment operator, and
   !! needs one coupled-perturbed solve per field direction.
   !!
   !! **The field is an imaginary, antisymmetric perturbation.** The
   !! first-order density it induces is antisymmetric, so it has no Coulomb
   !! response and, in Hartree-Fock, only exchange couples the equations.
   !! Everything here is the real matrix multiplying `i`.
   !!
   !! The working equations and their signs are `pyscf.prop.nmr.rhf`'s, and
   !! the tests are held to it.
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use mqc_physical_constants, only: FINE_STRUCTURE
   use mqc_calculation_defaults, only: PROPERTY_RESPONSE_TOL
   use mqc_scf_common, only: build_density_closed_shell
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_charges, only: ao_to_atom
   use mqc_math_utils, only: cross_product
   use mqc_czt_ao, only: eval_ao_block, AO_POINT_BLOCK
   use mqc_czt_xc, only: xc_context_t, xc_grid_gga_quantities
   use mqc_czt_pcm, only: pcm_context_t
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_response_product, only: response_mean_field
   use mqc_nmr_summary, only: shielding_summary
   use mqc_czt_giao, only: giao_one_electron, giao_coulomb_exchange, &
                           giao_smeared_charges, &
                           GIAO_OVLP, GIAO_KIN, GIAO_NUC, GIAO_RXP, GIAO_A11, &
                           GIAO_A01GP, GIAO_PRINVXP
   implicit none
   private

   public :: czt_nmr_shielding
   public :: print_shielding_report

   real(dp), parameter :: PER_MILLION = 1.0e6_dp
   real(dp), parameter :: TO_PPM = (FINE_STRUCTURE**2)*PER_MILLION
      !! A shielding in atomic units is in units of the fine-structure
      !! constant squared
   integer, parameter :: CP_MAX_ITER = 60
      !! Largest Krylov subspace the field response may build per direction
   integer, parameter :: LINE_LEN = 160
      !! Length of a formatted report line
   integer, parameter :: LABEL_LEN = 16
      !! Length of an atom's label in the report

contains

   subroutine czt_nmr_shielding(mol, coeff, orbital_energies, n_occ, shielding, error, &
                                diamagnetic, paramagnetic, xc, pcm)
      !! Shielding tensors of every nucleus, in parts per million
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)             !! (n_ao, n_mo), canonical
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ                    !! Doubly occupied orbitals
      real(dp), allocatable, intent(out) :: shielding(:, :, :)
         !! (3, 3, n_atoms). First index the field's direction, second the
         !! nuclear moment's. Not symmetric in general.
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: diamagnetic(:, :, :)
      real(dp), allocatable, intent(out), optional :: paramagnetic(:, :, :)
         !! The two parts, same shape and units. Only their sum is a physical
         !! quantity: the split depends on the gauge choice.
      type(xc_context_t), intent(inout), optional :: xc
         !! The reference's functional. Absent or inactive is Hartree-Fock.
         !! LDA, GGA and their global hybrids.
      type(pcm_context_t), intent(inout), optional :: pcm
         !! The continuum the reference was converged in. Absent or not
         !! enabled is the gas phase.

      real(dp), allocatable :: density(:, :), h10(:, :, :), s10(:, :, :), work(:, :, :)
      real(dp), allocatable :: vj(:, :, :), vk(:, :, :), bounds(:, :)
      real(dp), allocatable :: h1(:, :, :), s1(:, :, :), mo1(:, :, :), dm10(:, :, :)
      real(dp), allocatable :: dia(:, :, :), para(:, :, :), a11(:, :, :), a01(:, :, :)
      real(dp) :: tensor(3, 3)
      real(dp) :: trace, exchange
      integer :: n_ao, n_mo, n_vir, x, y, iatom, i
      logical :: kohn_sham

      if (error%has_error()) return
      n_ao = mol%nao
      n_mo = size(coeff, 2)
      n_vir = n_mo - n_occ
      if (n_occ < 1 .or. n_vir < 1) then
         call error%set(ERROR_VALIDATION, "a shielding needs both occupied and "// &
                        "virtual orbitals")
         return
      end if
      if (mol%necpbas > 0) then
         call error%set(ERROR_VALIDATION, "shieldings are not available with an "// &
                        "effective core potential: its field derivative is not built")
         return
      end if

      ! Hartree-Fock keeps all of its exchange and a functional the fraction
      ! it declares, in the first-order Fock matrix and in the response alike.
      kohn_sham = .false.
      if (present(xc)) kohn_sham = xc%active
      exchange = 1.0_dp
      if (kohn_sham) then
         exchange = xc%exx_fraction
         if (xc%any_mgga) then
            call error%set(ERROR_VALIDATION, "shieldings are not available with a "// &
                           "meta-GGA: its kinetic energy density depends on the "// &
                           "gauge, and that dependence is not built")
            return
         end if
         if (xc%range_separated) then
            call error%set(ERROR_VALIDATION, "shieldings are not available with a "// &
                           "range-separated functional: the field derivative of the "// &
                           "attenuated exchange integrals is not built")
            return
         end if
         if (abs(xc%nlc_b) > 0.0_dp .or. abs(xc%nlc_c) > 0.0_dp) then
            call error%set(ERROR_VALIDATION, "shieldings are not available with a "// &
                           "non-local correlation functional: its potential over "// &
                           "gauge-including orbitals is not built")
            return
         end if
      end if

      allocate (density(n_ao, n_ao))
      call build_density_closed_shell(coeff, n_occ, density)

      ! ---- the first-order Fock and overlap matrices -----------------------
      !
      ! h10 = -(1/2) <r_j x p> - J' + (1/2)(K' - K'^T) - <g nuc> - <g kin>,
      ! with J' and K' the two-electron integrals differentiated through the
      ! gauge factor, and s10 = -<g | >.
      call giao_one_electron(mol, GIAO_RXP, h10, error)
      if (error%has_error()) return
      h10 = -0.5_dp*h10
      call giao_coulomb_exchange(mol, density, vj, vk, error)
      if (error%has_error()) return
      do x = 1, 3
         h10(:, :, x) = h10(:, :, x) - vj(:, :, x) &
                        + 0.5_dp*exchange*(vk(:, :, x) - transpose(vk(:, :, x)))
      end do
      deallocate (vj, vk)
      if (kohn_sham) then
         call giao_xc_potential(xc, mol, density, work, error)
         if (error%has_error()) return
         h10 = h10 - work
      end if
      ! A continuum's surface charges answer the charge density, which the
      ! field does not change at first order: they stay as the reference
      ! left them, and their potential enters as the nuclei's does, through
      ! the gauge factor of the functions it is taken between.
      if (present(pcm)) then
         if (pcm%enabled) then
            block
               real(dp), allocatable :: q(:)

               call pcm%surface_charges(mol, density, q, error)
               if (error%has_error()) return
               call giao_smeared_charges(mol, pcm%surface%coords, pcm%surface%zeta, q, &
                                         work, error)
               if (error%has_error()) return
               h10 = h10 - work
            end block
         end if
      end if
      call giao_one_electron(mol, GIAO_NUC, work, error)
      if (error%has_error()) return
      h10 = h10 - work
      call giao_one_electron(mol, GIAO_KIN, work, error)
      if (error%has_error()) return
      h10 = h10 - work
      call giao_one_electron(mol, GIAO_OVLP, s10, error)
      if (error%has_error()) return
      s10 = -s10

      ! Both into the orbital basis, every orbital by the occupied ones.
      call occupied_columns(h10, coeff, n_occ, h1)
      call occupied_columns(s10, coeff, n_occ, s1)
      deallocate (h10, s10)

      ! ---- the orbitals' response to the field -----------------------------
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      call solve_field_response(mol, coeff, orbital_energies, n_occ, bounds, h1, s1, &
                                exchange, mo1, error)
      if (error%has_error()) return

      ! The first-order density's generator, twice for the two electrons an
      ! orbital holds: C (2 mo1) C_occ^T.
      call generator_density(mo1, coeff, n_occ, dm10)

      ! ---- the two parts, nucleus by nucleus -------------------------------
      allocate (dia(3, 3, mol%natm), para(3, 3, mol%natm))
      do iatom = 1, mol%natm
         call giao_one_electron(mol, GIAO_A11, a11, error, origin=mol%coords(:, iatom))
         call giao_one_electron(mol, GIAO_A01GP, a01, error, origin=mol%coords(:, iatom))
         call giao_one_electron(mol, GIAO_PRINVXP, work, error, origin=mol%coords(:, iatom))
         if (error%has_error()) return

         ! Diamagnetic: the nine components are a 3x3 tensor in row-major
         ! order, and the first integral enters with its trace removed.
         do x = 1, 3
            do y = 1, 3
               tensor(x, y) = sum(a11(:, :, 3*(x - 1) + y)*density)
            end do
         end do
         trace = tensor(1, 1) + tensor(2, 2) + tensor(3, 3)
         do i = 1, 3
            tensor(i, i) = tensor(i, i) - trace
         end do
         do x = 1, 3
            do y = 1, 3
               tensor(x, y) = tensor(x, y) + sum(a01(:, :, 3*(x - 1) + y)*density)
            end do
         end do
         dia(:, :, iatom) = tensor

         ! Paramagnetic: twice the first-order density against the
         ! nuclear-moment operator, the two for the complex conjugate.
         do x = 1, 3
            do y = 1, 3
               para(x, y, iatom) = 2.0_dp*sum(transpose(dm10(:, :, x))*work(:, :, y))
            end do
         end do
      end do

      dia = TO_PPM*dia
      para = TO_PPM*para
      shielding = dia + para
      if (present(diamagnetic)) call move_alloc(dia, diamagnetic)
      if (present(paramagnetic)) call move_alloc(para, paramagnetic)
   end subroutine czt_nmr_shielding

   subroutine solve_field_response(mol, coeff, eps, n_occ, bounds, h1, s1, exchange, &
                                   mo1, error)
      !! The imaginary first-order orbital coefficients, one set per direction
      !!
      !! The occupied-occupied block is fixed by the first-order overlap,
      !! `-s1/2`. The virtual-occupied block `u` solves
      !!
      !!     (e_a - e_i) u_ai + V_ai[u] = -(h1 - s1 e_i)_ai - V_ai[occupied block]
      !!
      !! where `V` is the exchange response to the antisymmetric density a
      !! set of coefficients generates. Solved in a Krylov subspace per
      !! direction, the three sharing each integral pass.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :), eps(:)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: bounds(:, :)
      real(dp), intent(in) :: h1(:, :, :), s1(:, :, :)   !! (n_mo, n_occ, 3)
      real(dp), intent(in) :: exchange
         !! The fraction of exact exchange the reference keeps. At zero, a
         !! pure functional, the field induces no potential -- a functional
         !! of the density alone does not see the current it drives -- and
         !! the equations are solved by the orbital-energy gaps in one step.
      real(dp), allocatable, intent(out) :: mo1(:, :, :)   !! (n_mo, n_occ, 3)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gaps(:), rhs(:, :), resid(:, :), trial(:, :, :), image(:, :, :)
      real(dp), allocatable :: small(:, :, :), proj(:, :), y(:), full(:, :, :), v(:, :, :)
      real(dp) :: rnorm(3)
      real(dp) :: norm
      integer :: n_done(3)
      integer :: n_mo, n_vir, n, x, i, a, k, m, pass
      logical :: active(3)
      logical :: singular
      character(len=LINE_LEN) :: line

      n_mo = size(coeff, 2)
      n_vir = n_mo - n_occ
      n = n_vir*n_occ
      allocate (mo1(n_mo, n_occ, 3))
      mo1 = 0.0_dp
      mo1(1:n_occ, :, :) = -0.5_dp*s1(1:n_occ, :, :)

      allocate (gaps(n))
      do i = 1, n_occ
         do a = 1, n_vir
            gaps(a + (i - 1)*n_vir) = eps(n_occ + a) - eps(i)
         end do
      end do

      ! Right-hand side: the perturbation, and the response to the fixed
      ! occupied-occupied block.
      call exchange_response(mol, coeff, n_occ, bounds, exchange, mo1, v, error)
      if (error%has_error()) return
      allocate (rhs(n, 3), resid(n, 3))
      do x = 1, 3
         do i = 1, n_occ
            do a = 1, n_vir
               rhs(a + (i - 1)*n_vir, x) = -(h1(n_occ + a, i, x) - s1(n_occ + a, i, x)*eps(i)) &
                                           - v(n_occ + a, i, x)
            end do
         end do
      end do

      allocate (trial(n, CP_MAX_ITER, 3), image(n, CP_MAX_ITER, 3))
      allocate (small(CP_MAX_ITER, CP_MAX_ITER, 3), proj(CP_MAX_ITER, 3))
      allocate (full(n_mo, n_occ, 3))
      resid = rhs
      active = .true.
      n_done = 0
      do k = 1, min(CP_MAX_ITER, n)
         do x = 1, 3
            rnorm(x) = largest_residual(resid(:, x))
            if (rnorm(x) < PROPERTY_RESPONSE_TOL) active(x) = .false.
         end do
         if (.not. any(active)) exit

         full = 0.0_dp
         do x = 1, 3
            if (.not. active(x)) cycle
            trial(:, k, x) = resid(:, x)/gaps
            do pass = 1, 2
               do m = 1, k - 1
                  trial(:, k, x) = trial(:, k, x) &
                                   - dot_product(trial(:, m, x), trial(:, k, x))*trial(:, m, x)
               end do
            end do
            norm = sqrt(dot_product(trial(:, k, x), trial(:, k, x)))
            if (norm < tiny(1.0_dp)) then
               active(x) = .false.
               cycle
            end if
            trial(:, k, x) = trial(:, k, x)/norm
            full(n_occ + 1:, :, x) = reshape(trial(:, k, x), [n_vir, n_occ])
         end do
         if (.not. any(active)) exit

         call exchange_response(mol, coeff, n_occ, bounds, exchange, full, v, error)
         if (error%has_error()) return
         do x = 1, 3
            if (.not. active(x)) cycle
            image(:, k, x) = gaps*trial(:, k, x) + reshape(v(n_occ + 1:, :, x), [n])
            do m = 1, k
               small(m, k, x) = dot_product(trial(:, m, x), image(:, k, x))
               small(k, m, x) = dot_product(trial(:, k, x), image(:, m, x))
            end do
            proj(k, x) = dot_product(trial(:, k, x), rhs(:, x))
            call solve_small(small(1:k, 1:k, x), proj(1:k, x), y, singular)
            if (singular) then
               call error%set(ERROR_VALIDATION, "the shielding's coupled-perturbed "// &
                              "equations are singular")
               return
            end if
            resid(:, x) = rhs(:, x) - matmul(image(:, 1:k, x), y)
            mo1(n_occ + 1:, :, x) = reshape(matmul(trial(:, 1:k, x), y), [n_vir, n_occ])
            n_done(x) = k
         end do
      end do
      do x = 1, 3
         rnorm(x) = largest_residual(resid(:, x))
      end do
      ! Written so that a residual that is not a number fails the test.
      if (.not. all(rnorm < PROPERTY_RESPONSE_TOL)) then
         call error%set(ERROR_VALIDATION, "the shielding's coupled-perturbed equations "// &
                        "did not converge")
         return
      end if
      write (line, "(a,i0,a,es9.2)") "  field response: ", maxval(n_done), &
         " iterations, residual ", maxval(rnorm)
      call logger%verbose(trim(line))
   end subroutine solve_field_response

   subroutine exchange_response(mol, coeff, n_occ, bounds, exchange, mo1, v, error)
      !! The exchange response to a set of first-order coefficients
      !!
      !! The density is `M - M^T` with `M = C (2 mo1) C_occ^T`: antisymmetric,
      !! so it has no Coulomb term and the mean field is `-K/2` alone. What
      !! comes back is that matrix between every orbital and the occupied
      !! ones.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: bounds(:, :)
      real(dp), intent(in) :: exchange   !! Fraction of exact exchange; zero returns zero
      real(dp), intent(in) :: mo1(:, :, :)                 !! (n_mo, n_occ, 3)
      real(dp), allocatable, intent(out) :: v(:, :, :)     !! (n_mo, n_occ, 3)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dens(:, :, :), g(:, :, :), zero_h(:, :)
      integer :: n_ao, x

      n_ao = size(coeff, 1)
      if (abs(exchange) <= 0.0_dp) then
         allocate (v(size(coeff, 2), n_occ, 3))
         v = 0.0_dp
         return
      end if
      allocate (zero_h(n_ao, n_ao))
      zero_h = 0.0_dp
      call generator_density(mo1, coeff, n_occ, dens)
      do x = 1, 3
         dens(:, :, x) = dens(:, :, x) - transpose(dens(:, :, x))
      end do
      call response_mean_field(mol, dens, zero_h, g, error, minus=.true., bounds=bounds, &
                               k_scale=exchange)
      if (error%has_error()) return
      call occupied_columns(g, coeff, n_occ, v)
   end subroutine exchange_response

   subroutine giao_xc_potential(xc, mol, density, vmat, error)
      !! The exchange-correlation potential's first-order matrix in the field
      !!
      !! The potential itself does not respond: it is a functional of the
      !! density, which a magnetic field leaves alone at first order. What
      !! is left is the gauge factor each basis function carries, whose field
      !! derivative at a point `r` for a function centred at `R` is
      !! `(1/2) (R x r) chi(r)`, in units of `i`. So the matrix is the
      !! potential between a function and another's derivative, less its
      !! transpose, and for a GGA the gradient of that derivative as well.
      !! `get_vxc_giao` of `pyscf.prop.nmr.rks`.
      type(xc_context_t), intent(inout) :: xc
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! Total, closed shell
      real(dp), allocatable, intent(out) :: vmat(:, :, :)
         !! (n_ao, n_ao, 3), antisymmetric, last index the field's direction
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), exc(:), vrho(:), gcoef(:, :)
      real(dp), allocatable :: ao(:, :), aog(:, :, :), aow(:, :), ig(:, :), dig(:, :)
      real(dp), allocatable :: part(:, :), axes(:, :)
      integer, allocatable :: owner(:)
      real(dp) :: r(3), cross(3), wv(3)
      real(dp) :: w0
      integer :: n_ao, npts, g0, g1, nb, p, u, a, b

      if (error%has_error()) return
      n_ao = mol%nao
      npts = xc%grid%n_points
      call xc_grid_gga_quantities(xc, mol, density, rho, exc, vrho, gcoef, error)
      if (error%has_error()) return

      ! Which atom each basis function sits on.
      call ao_to_atom(mol, owner)

      allocate (vmat(n_ao, n_ao, 3), part(n_ao, n_ao), axes(3, n_ao))
      vmat = 0.0_dp
      do g0 = 1, npts, AO_POINT_BLOCK
         g1 = min(g0 + AO_POINT_BLOCK - 1, npts)
         nb = g1 - g0 + 1
         if (xc%any_gga) then
            call eval_ao_block(mol, xc%grid%coords(:, g0:g1), ao, error, grad=aog)
         else
            call eval_ao_block(mol, xc%grid%coords(:, g0:g1), ao, error)
         end if
         if (error%has_error()) return
         if (allocated(aow)) then
            if (size(aow, 1) /= nb) deallocate (aow, ig, dig)
         end if
         if (.not. allocated(aow)) allocate (aow(nb, n_ao), ig(nb, n_ao), dig(nb, n_ao))

         ! The function weighted by the potential: `w v chi`, and for a GGA
         ! `w (dE/d grad rho) . grad chi` beside it.
         do u = 1, n_ao
            do p = 1, nb
               aow(p, u) = xc%grid%weights(g0 + p - 1)*vrho(g0 + p - 1)*ao(p, u)
            end do
         end do
         if (xc%any_gga) then
            do a = 1, 3
               do u = 1, n_ao
                  do p = 1, nb
                     aow(p, u) = aow(p, u) + xc%grid%weights(g0 + p - 1) &
                                 *gcoef(g0 + p - 1, a)*aog(p, u, a)
                  end do
               end do
            end do
         end if

         do b = 1, 3
            ! The gauge derivative of every function along this direction,
            ! and for a GGA its gradient contracted with `w dE/d grad rho`.
            ! d/dr_a of (R x r)_b is (R x e_a)_b, the same at every point.
            do u = 1, n_ao
               axes(:, u) = cross_with_axes(mol%coords(:, owner(u)), b)
            end do
            do u = 1, n_ao
               do p = 1, nb
                  r = xc%grid%coords(:, g0 + p - 1)
                  cross = cross_product(mol%coords(:, owner(u)), r)
                  ig(p, u) = 0.5_dp*cross(b)*ao(p, u)
                  if (xc%any_gga) then
                     w0 = xc%grid%weights(g0 + p - 1)
                     wv = w0*gcoef(g0 + p - 1, :)
                     dig(p, u) = 0.5_dp*(dot_product(wv, axes(:, u)) &
                                         *ao(p, u) &
                                         + cross(b)*(wv(1)*aog(p, u, 1) + wv(2)*aog(p, u, 2) &
                                                     + wv(3)*aog(p, u, 3)))
                  end if
               end do
            end do
            call pic_gemm(aow, ig, part, transa="T", beta=0.0_dp)
            vmat(:, :, b) = vmat(:, :, b) + part
            if (xc%any_gga) then
               call pic_gemm(ao, dig, part, transa="T", beta=0.0_dp)
               vmat(:, :, b) = vmat(:, :, b) + part
            end if
         end do
      end do
      do b = 1, 3
         vmat(:, :, b) = vmat(:, :, b) - transpose(vmat(:, :, b))
      end do
   end subroutine giao_xc_potential

   pure function cross_with_axes(a, component) result(c)
      !! `(a x e_k)_component` for the three unit vectors `e_k`
      real(dp), intent(in) :: a(3)
      integer, intent(in) :: component
      real(dp) :: c(3)

      real(dp) :: unit(3), full(3)
      integer :: k

      do k = 1, 3
         unit = 0.0_dp
         unit(k) = 1.0_dp
         full = cross_product(a, unit)
         c(k) = full(component)
      end do
   end function cross_with_axes

   subroutine occupied_columns(ao, coeff, n_occ, mo)
      !! `C^T A C_occ` for each of a stack of matrices in the atomic basis
      real(dp), intent(in) :: ao(:, :, :)        !! (n_ao, n_ao, n)
      real(dp), intent(in) :: coeff(:, :)        !! (n_ao, n_mo)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: mo(:, :, :)   !! (n_mo, n_occ, n)

      real(dp), allocatable :: half(:, :)
      integer :: x

      allocate (mo(size(coeff, 2), n_occ, size(ao, 3)), half(size(coeff, 1), n_occ))
      do x = 1, size(ao, 3)
         call pic_gemm(ao(:, :, x), coeff(:, 1:n_occ), half, beta=0.0_dp)
         call pic_gemm(coeff, half, mo(:, :, x), transa="T", beta=0.0_dp)
      end do
   end subroutine occupied_columns

   subroutine generator_density(mo1, coeff, n_occ, dens)
      !! `C (2 mo1) C_occ^T`: what a set of first-order coefficients adds to
      !! the density, before it is antisymmetrised
      real(dp), intent(in) :: mo1(:, :, :)       !! (n_mo, n_occ, n)
      real(dp), intent(in) :: coeff(:, :)        !! (n_ao, n_mo)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: dens(:, :, :)   !! (n_ao, n_ao, n)

      real(dp), allocatable :: half(:, :)
      integer :: n_ao, x

      n_ao = size(coeff, 1)
      allocate (dens(n_ao, n_ao, size(mo1, 3)), half(n_ao, n_occ))
      do x = 1, size(mo1, 3)
         call pic_gemm(coeff, mo1(:, :, x), half, alpha=2.0_dp, beta=0.0_dp)
         call pic_gemm(half, coeff(:, 1:n_occ), dens(:, :, x), transb="T", beta=0.0_dp)
      end do
   end subroutine generator_density

   pure function largest_residual(resid) result(norm)
      !! Largest element of a residual in magnitude; `huge` when any element
      !! is not a number, so that it can never read as converged
      real(dp), intent(in) :: resid(:)
      real(dp) :: norm

      if (ieee_is_nan(sum(resid))) then
         norm = huge(1.0_dp)
      else
         norm = maxval(abs(resid))
      end if
   end function largest_residual

   subroutine solve_small(matrix, rhs, x, singular)
      !! A small dense system by LU factorisation
      real(dp), intent(in) :: matrix(:, :), rhs(:)
      real(dp), allocatable, intent(out) :: x(:)
      logical, intent(out) :: singular   !! The factorisation met a zero pivot

      real(dp), allocatable :: a(:, :), b(:, :)
      integer, allocatable :: ipiv(:)
      integer :: n, info

      n = size(rhs)
      allocate (a(n, n), source=matrix)
      allocate (b(n, 1), ipiv(n))
      b(:, 1) = rhs
      call pic_getrf(a, ipiv, info)
      singular = info /= 0
      if (singular) return
      call pic_getrs(a, ipiv, b, info=info)
      singular = info /= 0
      allocate (x(n), source=b(:, 1))
   end subroutine solve_small

   subroutine print_shielding_report(shielding, element_symbols)
      !! Isotropic shielding and anisotropy per nucleus, then the tensors
      real(dp), intent(in) :: shielding(:, :, :)   !! (3, 3, n_atoms), ppm
      character(len=*), intent(in) :: element_symbols(:)

      character(len=LINE_LEN) :: line
      character(len=LABEL_LEN) :: label
      real(dp) :: iso, aniso
      integer :: iatom, x

      call logger%info("")
      call logger%info("  ==== nuclear magnetic shielding (GIAO), ppm ====")
      call logger%info("")
      call logger%info("     atom       isotropic    anisotropy")
      do iatom = 1, size(shielding, 3)
         call shielding_summary(shielding(:, :, iatom), iso, aniso)
         write (label, "(a,i0)") trim(adjustl(element_symbols(iatom)))//" ", iatom
         write (line, "(4x,a8,2f14.4)") label, iso, aniso
         call logger%info(trim(line))
      end do
      call logger%verbose("")
      do iatom = 1, size(shielding, 3)
         write (label, "(a,i0)") trim(adjustl(element_symbols(iatom)))//" ", iatom
         call logger%verbose("    "//trim(label)//", field direction by row")
         do x = 1, 3
            write (line, "(6x,3f14.4)") shielding(x, :, iatom)
            call logger%verbose(trim(line))
         end do
      end do
      call logger%info("")
   end subroutine print_shielding_report

end module mqc_czt_nmr
