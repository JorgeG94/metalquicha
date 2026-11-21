!! Integrals over gauge-including atomic orbitals, for magnetic properties
module mqc_czt_giao
   !! The matrices a nuclear magnetic shielding needs, each the first
   !! derivative of an ordinary integral with respect to an external magnetic
   !! field when every basis function carries its own gauge origin (a London
   !! orbital), or with respect to a nuclear magnetic moment.
   !!
   !! All are imaginary in the complex basis the field induces; what is
   !! returned is the real matrix multiplying `i`, under libcint's sign for
   !! each integral, which is what the formulas in `mqc_czt_nmr` are written
   !! against.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim
   use libcint_fortran, only: LIBCINT_PTR_RINV_ORIG, &
                              libcint_1e_igovlp_sph, libcint_1e_igovlp_cart, &
                              libcint_1e_igkin_sph, libcint_1e_igkin_cart, &
                              libcint_1e_ignuc_sph, libcint_1e_ignuc_cart, &
                              libcint_1e_giao_irjxp_sph, libcint_1e_giao_irjxp_cart, &
                              libcint_1e_giao_a11part_sph, libcint_1e_giao_a11part_cart, &
                              libcint_1e_a01gp_sph, libcint_1e_a01gp_cart, &
                              libcint_1e_prinvxp_sph, libcint_1e_prinvxp_cart, &
                              libcint_2e_ig1_sph, libcint_2e_ig1_cart
   implicit none
   private

   public :: giao_one_electron
   public :: giao_coulomb_exchange
   public :: GIAO_OVLP, GIAO_KIN, GIAO_NUC, GIAO_RXP, GIAO_A11, GIAO_A01GP, GIAO_PRINVXP

   integer, parameter :: GIAO_OVLP = 1     !! `int1e_igovlp`, 3 components
   integer, parameter :: GIAO_KIN = 2      !! `int1e_igkin`, 3
   integer, parameter :: GIAO_NUC = 3      !! `int1e_ignuc`, 3
   integer, parameter :: GIAO_RXP = 4      !! `int1e_giao_irjxp`, 3
   integer, parameter :: GIAO_A11 = 5      !! `int1e_giao_a11part`, 9, needs an origin
   integer, parameter :: GIAO_A01GP = 6    !! `int1e_a01gp`, 9, needs an origin
   integer, parameter :: GIAO_PRINVXP = 7  !! `int1e_prinvxp`, 3, needs an origin
   integer, parameter :: N_VECTOR = 3
      !! Components of an integral that is a vector in the field
   integer, parameter :: N_TENSOR = 9
      !! Components of one that is a tensor in field and nuclear moment

contains

   pure function component_count(which) result(n)
      !! How many components the integral `which` has
      integer, intent(in) :: which   !! One of the `GIAO_*` constants
      integer :: n

      select case (which)
      case (GIAO_A11, GIAO_A01GP)
         n = N_TENSOR
      case default
         n = N_VECTOR
      end select
   end function component_count

   function max_shell_dim(mol) result(mx)
      !! The largest number of functions in one shell of `mol`
      type(czt_molecule_t), intent(in) :: mol
      integer :: mx

      integer :: ish

      mx = 0
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do
   end function max_shell_dim

   subroutine giao_one_electron(mol, which, matrix, error, origin)
      !! One of the one-electron matrices, every ordered shell pair computed
      !!
      !! No symmetry is assumed between `(i, j)` and `(j, i)`: some of these
      !! are antisymmetric, some neither, and the full square costs little.
      type(czt_molecule_t), intent(in) :: mol
      integer, intent(in) :: which   !! One of the `GIAO_*` constants
      real(dp), allocatable, intent(out) :: matrix(:, :, :)
         !! (n_ao, n_ao, n_comp). The nine-component ones are a 3x3 tensor in
         !! row-major order: component `3*(x-1) + y`.
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: origin(3)
         !! Where the `1/r` operator sits, in Bohr: the nucleus whose
         !! shielding is wanted. Required by the three that carry one.

      real(dp), allocatable :: env(:), buf(:)
      integer :: shls(2)
      integer :: n_comp, ish, jsh, di, dj, io, jo, i, j, comp, mx
      logical :: have

      if (error%has_error()) return
      n_comp = component_count(which)
      select case (which)
      case (GIAO_A11, GIAO_A01GP, GIAO_PRINVXP)
         if (.not. present(origin)) then
            call error%set(ERROR_VALIDATION, "a nuclear-moment integral was asked "// &
                           "for without the nucleus it belongs to")
            return
         end if
      case (GIAO_OVLP, GIAO_KIN, GIAO_NUC, GIAO_RXP)
         continue
      case default
         call error%set(ERROR_VALIDATION, "unknown gauge-including integral")
         return
      end select

      allocate (env(size(mol%env)), source=mol%env)
      if (present(origin)) env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = origin

      allocate (matrix(mol%nao, mol%nao, n_comp))
      matrix = 0.0_dp
      mx = max_shell_dim(mol)

      !$omp parallel default(none) &
      !$omp    shared(mol, matrix, env, mx, n_comp, which) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, shls, buf, have)
      allocate (buf(mx*mx*n_comp))
      !$omp do schedule(dynamic)
      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]
            have = one_shell_pair(which, mol, buf, shls, env)
            if (.not. have) cycle
            do comp = 1, n_comp
               do j = 1, dj
                  do i = 1, di
                     matrix(io + i, jo + j, comp) = buf(i + di*(j - 1 + dj*(comp - 1)))
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf)
      !$omp end parallel
   end subroutine giao_one_electron

   function one_shell_pair(which, mol, buf, shls, env) result(have)
      !! One shell pair of the integral `which`; false when it vanishes
      integer, intent(in) :: which   !! One of the `GIAO_*` constants
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(inout) :: buf(:)
         !! (d_i, d_j, n_comp) flattened, first index fastest
      integer, intent(in) :: shls(2)   !! The two shells, 0-based
      real(dp), intent(in) :: env(:)
         !! `mol%env`, or a copy with the `1/r` origin set
      logical :: have

      integer :: ret

      if (mol%cartesian) then
         select case (which)
         case (GIAO_OVLP)
            ret = libcint_1e_igovlp_cart(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_KIN)
            ret = libcint_1e_igkin_cart(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_NUC)
            ret = libcint_1e_ignuc_cart(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_RXP)
            ret = libcint_1e_giao_irjxp_cart(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_A11)
            ret = libcint_1e_giao_a11part_cart(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_A01GP)
            ret = libcint_1e_a01gp_cart(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_PRINVXP)
            ret = libcint_1e_prinvxp_cart(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case default
            ret = 0
         end select
      else
         select case (which)
         case (GIAO_OVLP)
            ret = libcint_1e_igovlp_sph(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_KIN)
            ret = libcint_1e_igkin_sph(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_NUC)
            ret = libcint_1e_ignuc_sph(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_RXP)
            ret = libcint_1e_giao_irjxp_sph(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_A11)
            ret = libcint_1e_giao_a11part_sph(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_A01GP)
            ret = libcint_1e_a01gp_sph(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case (GIAO_PRINVXP)
            ret = libcint_1e_prinvxp_sph(buf, shls, mol%atm, mol%natm, mol%bas, mol%nbas, env)
         case default
            ret = 0
         end select
      end if
      have = ret /= 0
   end function one_shell_pair

   subroutine giao_coulomb_exchange(mol, density, vj, vk, error)
      !! The field derivative of the two-electron integrals against a density
      !!
      !! With `(g i j|k l)` the integral differentiated through the gauge
      !! factor of the first pair,
      !!
      !!     vj_ij = sum_kl (g i j|k l) D_lk,      vk_il = sum_jk (g i j|k l) D_jk
      !!
      !! for each of the three field directions. Neither is symmetrised or
      !! signed here; the caller applies the convention it wants.
      ! TODO(mqc): every shell quartet is computed, with no Schwarz screening
      ! and no use of the symmetry of the second pair, so this is the full
      ! fourth power of the basis size and the slowest step of a shielding.
      ! A bound on the differentiated integral would let it screen as the
      ! Fock build does, and `(kl) = (lk)` would halve what is left.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! Symmetric
      real(dp), allocatable, intent(out) :: vj(:, :, :), vk(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: buf(:), vj_t(:, :, :), vk_t(:, :, :)
      integer :: shls(4)
      integer :: ret
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo
      integer :: i, j, k, l, comp, mx, idx, n
      real(dp) :: val

      if (error%has_error()) return
      n = mol%nao
      mx = max_shell_dim(mol)
      allocate (vj(n, n, N_VECTOR), vk(n, n, N_VECTOR))
      vj = 0.0_dp
      vk = 0.0_dp

      !$omp parallel default(none) shared(mol, density, vj, vk, mx, n) &
      !$omp    private(ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, i, j, k, l, &
      !$omp            comp, idx, val, shls, buf, ret, vj_t, vk_t)
      allocate (buf(N_VECTOR*mx**4), vj_t(n, n, N_VECTOR), vk_t(n, n, N_VECTOR))
      vj_t = 0.0_dp
      vk_t = 0.0_dp
      !$omp do schedule(dynamic) collapse(2)
      do ish = 1, mol%nbas
         do jsh = 1, mol%nbas
            di = shell_dim(mol%cartesian, ish - 1, mol%bas)
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            io = mol%shell_offset(ish)
            jo = mol%shell_offset(jsh)
            do ksh = 1, mol%nbas
               dk = shell_dim(mol%cartesian, ksh - 1, mol%bas)
               ko = mol%shell_offset(ksh)
               do lsh = 1, mol%nbas
                  dl = shell_dim(mol%cartesian, lsh - 1, mol%bas)
                  lo = mol%shell_offset(lsh)
                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]
                  if (mol%cartesian) then
                     ret = libcint_2e_ig1_cart(buf, shls, mol%atm, mol%natm, &
                                               mol%bas, mol%nbas, mol%env)
                  else
                     ret = libcint_2e_ig1_sph(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, mol%env)
                  end if
                  if (ret == 0) cycle
                  idx = 0
                  do comp = 1, N_VECTOR
                     do l = 1, dl
                        do k = 1, dk
                           do j = 1, dj
                              do i = 1, di
                                 idx = idx + 1
                                 val = buf(idx)
                                 vj_t(io + i, jo + j, comp) = vj_t(io + i, jo + j, comp) &
                                                              + val*density(lo + l, ko + k)
                                 vk_t(io + i, lo + l, comp) = vk_t(io + i, lo + l, comp) &
                                                              + val*density(jo + j, ko + k)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      !$omp critical
      vj = vj + vj_t
      vk = vk + vk_t
      !$omp end critical
      deallocate (buf, vj_t, vk_t)
      !$omp end parallel
   end subroutine giao_coulomb_exchange

end module mqc_czt_giao
