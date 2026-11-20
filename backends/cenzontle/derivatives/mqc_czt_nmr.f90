!! Nuclear magnetic shielding tensors for a closed-shell Hartree-Fock reference
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
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_response_product, only: response_mean_field
   use mqc_nmr_summary, only: shielding_summary
   use mqc_czt_giao, only: giao_one_electron, giao_coulomb_exchange, &
                           GIAO_OVLP, GIAO_KIN, GIAO_NUC, GIAO_RXP, GIAO_A11, &
                           GIAO_A01GP, GIAO_PRINVXP
   implicit none
   private

   public :: czt_nmr_shielding
   public :: print_shielding_report
   public :: FINE_STRUCTURE

   real(dp), parameter :: FINE_STRUCTURE = 7.2973525693e-3_dp
      !! The fine-structure constant, CODATA 2018. A shielding in atomic units
      !! is in units of its square.
   real(dp), parameter :: TO_PPM = FINE_STRUCTURE**2*1.0e6_dp
   integer, parameter :: CP_MAX_ITER = 60
   real(dp), parameter :: CP_TOL = 1.0e-10_dp
      !! Largest element of the coupled-perturbed residual at convergence
   real(dp), parameter :: SINGULAR_PIVOT = 1.0e-13_dp
   integer, parameter :: LINE_LEN = 160

contains

   subroutine czt_nmr_shielding(mol, coeff, orbital_energies, n_occ, shielding, error, &
                                diamagnetic, paramagnetic)
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

      real(dp), allocatable :: density(:, :), h10(:, :, :), s10(:, :, :), work(:, :, :)
      real(dp), allocatable :: vj(:, :, :), vk(:, :, :), bounds(:, :)
      real(dp), allocatable :: h1(:, :, :), s1(:, :, :), mo1(:, :, :), dm10(:, :, :)
      real(dp), allocatable :: dia(:, :, :), para(:, :, :), a11(:, :, :), a01(:, :, :)
      real(dp), allocatable :: half(:, :)
      real(dp) :: tensor(3, 3)
      real(dp) :: trace
      integer :: n_ao, n_mo, n_vir, x, y, iatom, i

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

      allocate (density(n_ao, n_ao))
      call pic_gemm(coeff(:, 1:n_occ), coeff(:, 1:n_occ), density, transb="T", &
                    alpha=2.0_dp, beta=0.0_dp)

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
                        + 0.5_dp*(vk(:, :, x) - transpose(vk(:, :, x)))
      end do
      deallocate (vj, vk)
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
      allocate (h1(n_mo, n_occ, 3), s1(n_mo, n_occ, 3), half(n_ao, n_occ))
      do x = 1, 3
         call pic_gemm(h10(:, :, x), coeff(:, 1:n_occ), half, beta=0.0_dp)
         call pic_gemm(coeff, half, h1(:, :, x), transa="T", beta=0.0_dp)
         call pic_gemm(s10(:, :, x), coeff(:, 1:n_occ), half, beta=0.0_dp)
         call pic_gemm(coeff, half, s1(:, :, x), transa="T", beta=0.0_dp)
      end do
      deallocate (h10, s10, half)

      ! ---- the orbitals' response to the field -----------------------------
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      call solve_field_response(mol, coeff, orbital_energies, n_occ, bounds, h1, s1, &
                                mo1, error)
      if (error%has_error()) return

      ! The first-order density's generator, twice for the two electrons an
      ! orbital holds: C (2 mo1) C_occ^T.
      allocate (dm10(n_ao, n_ao, 3), half(n_ao, n_occ))
      do x = 1, 3
         call pic_gemm(coeff, mo1(:, :, x), half, alpha=2.0_dp, beta=0.0_dp)
         call pic_gemm(half, coeff(:, 1:n_occ), dm10(:, :, x), transb="T", beta=0.0_dp)
      end do

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

   subroutine solve_field_response(mol, coeff, eps, n_occ, bounds, h1, s1, mo1, error)
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
      call exchange_response(mol, coeff, n_occ, bounds, mo1, v, error)
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
            rnorm(x) = maxval(abs(resid(:, x)))
            if (rnorm(x) < CP_TOL) active(x) = .false.
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

         call exchange_response(mol, coeff, n_occ, bounds, full, v, error)
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
         rnorm(x) = maxval(abs(resid(:, x)))
      end do
      if (maxval(rnorm) >= CP_TOL) then
         call error%set(ERROR_VALIDATION, "the shielding's coupled-perturbed equations "// &
                        "did not converge")
         return
      end if
      write (line, "(a,i0,a,es9.2)") "  field response: ", maxval(n_done), &
         " iterations, residual ", maxval(rnorm)
      call logger%verbose(trim(line))
   end subroutine solve_field_response

   subroutine exchange_response(mol, coeff, n_occ, bounds, mo1, v, error)
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
      real(dp), intent(in) :: mo1(:, :, :)                 !! (n_mo, n_occ, 3)
      real(dp), allocatable, intent(out) :: v(:, :, :)     !! (n_mo, n_occ, 3)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dens(:, :, :), g(:, :, :), half(:, :), zero_h(:, :)
      integer :: n_ao, n_mo, x

      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)
      allocate (dens(n_ao, n_ao, 3), half(n_ao, n_occ), zero_h(n_ao, n_ao))
      zero_h = 0.0_dp
      do x = 1, 3
         call pic_gemm(coeff, mo1(:, :, x), half, alpha=2.0_dp, beta=0.0_dp)
         call pic_gemm(half, coeff(:, 1:n_occ), dens(:, :, x), transb="T", beta=0.0_dp)
         dens(:, :, x) = dens(:, :, x) - transpose(dens(:, :, x))
      end do
      call response_mean_field(mol, dens, zero_h, g, error, minus=.true., bounds=bounds, &
                               k_scale=1.0_dp)
      if (error%has_error()) return
      allocate (v(n_mo, n_occ, 3))
      do x = 1, 3
         call pic_gemm(g(:, :, x), coeff(:, 1:n_occ), half, beta=0.0_dp)
         call pic_gemm(coeff, half, v(:, :, x), transa="T", beta=0.0_dp)
      end do
   end subroutine exchange_response

   subroutine solve_small(matrix, rhs, x, singular)
      !! A small dense system by elimination with partial pivoting
      real(dp), intent(in) :: matrix(:, :), rhs(:)
      real(dp), allocatable, intent(out) :: x(:)
      logical, intent(out) :: singular

      real(dp), allocatable :: a(:, :), row(:)
      real(dp) :: scale, factor
      integer :: n, i, j, pivot

      n = size(rhs)
      a = matrix
      x = rhs
      singular = .false.
      scale = maxval(abs(a))
      allocate (row(n))
      do i = 1, n
         pivot = i - 1 + maxloc(abs(a(i:n, i)), dim=1)
         if (abs(a(pivot, i)) <= SINGULAR_PIVOT*scale) then
            singular = .true.
            return
         end if
         if (pivot /= i) then
            row = a(i, :)
            a(i, :) = a(pivot, :)
            a(pivot, :) = row
            factor = x(i)
            x(i) = x(pivot)
            x(pivot) = factor
         end if
         do j = i + 1, n
            factor = a(j, i)/a(i, i)
            a(j, i:n) = a(j, i:n) - factor*a(i, i:n)
            x(j) = x(j) - factor*x(i)
         end do
      end do
      do i = n, 1, -1
         x(i) = (x(i) - dot_product(a(i, i + 1:n), x(i + 1:n)))/a(i, i)
      end do
   end subroutine solve_small

   subroutine print_shielding_report(shielding, element_symbols)
      !! Isotropic shielding and anisotropy per nucleus, then the tensors
      real(dp), intent(in) :: shielding(:, :, :)   !! (3, 3, n_atoms), ppm
      character(len=*), intent(in) :: element_symbols(:)

      character(len=LINE_LEN) :: line
      character(len=16) :: label
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
