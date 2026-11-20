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
   ! These integrals are generated in libfint and are not among the entry
   ! points its `libcint_fortran` module exports, so this module reaches the
   ! generated functions themselves. It is the only one that does: when
   ! `libcint_fortran` exports them, this file is what changes.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim
   use libcint_fortran, only: LIBCINT_PTR_RINV_ORIG
   use cint_workspace, only: cint_ws
   use cint_gen_intor1, only: int1e_igovlp_sph, int1e_igovlp_cart, &
                              int1e_igkin_sph, int1e_igkin_cart, &
                              int1e_ignuc_sph, int1e_ignuc_cart, &
                              int1e_giao_irjxp_sph, int1e_giao_irjxp_cart, &
                              int1e_giao_a11part_sph, int1e_giao_a11part_cart, &
                              int1e_a01gp_sph, int1e_a01gp_cart, &
                              int1e_prinvxp_sph, int1e_prinvxp_cart
   use cint_gen_intor2, only: int2e_ig1_sph, int2e_ig1_cart
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

contains

   pure function component_count(which) result(n)
      integer, intent(in) :: which
      integer :: n

      select case (which)
      case (GIAO_A11, GIAO_A01GP)
         n = 9
      case default
         n = 3
      end select
   end function component_count

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

      integer, allocatable :: atm(:), bas(:)
      real(dp), allocatable :: env(:), buf(:)
      type(cint_ws) :: ws
      integer :: dims(0:3), shls(0:3)
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

      allocate (atm(size(mol%atm)), bas(size(mol%bas)))
      atm = reshape(mol%atm, [size(mol%atm)])
      bas = reshape(mol%bas, [size(mol%bas)])
      allocate (env(size(mol%env)), source=mol%env)
      if (present(origin)) env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = origin

      allocate (matrix(mol%nao, mol%nao, n_comp))
      matrix = 0.0_dp
      mx = 0
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do

      !$omp parallel default(none) &
      !$omp    shared(mol, matrix, atm, bas, env, mx, n_comp, which) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, shls, dims, buf, ws, have)
      ! A private copy of the workspace is not default-initialised, and the
      ! integral driver reads this pointer before anything else.
      ws%opt => null()
      allocate (buf(mx*mx*n_comp))
      !$omp do schedule(dynamic)
      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1, 0, 0]
            dims = [di, dj, 1, 1]
            have = one_shell_pair(which, mol%cartesian, buf(1:di*dj*n_comp), dims, shls, &
                                  atm, mol%natm, bas, mol%nbas, env, ws)
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

   function one_shell_pair(which, cartesian, buf, dims, shls, atm, natm, bas, nbas, &
                           env, ws) result(have)
      !! Dispatch one shell pair to the generated integral
      integer, intent(in) :: which
      logical, intent(in) :: cartesian
      real(dp), intent(inout) :: buf(0:)
      integer, intent(in) :: dims(0:), shls(0:), atm(0:), bas(0:)
      integer, intent(in) :: natm, nbas
      real(dp), intent(in) :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: have

      have = .false.
      if (cartesian) then
         select case (which)
         case (GIAO_OVLP)
            have = int1e_igovlp_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_KIN)
            have = int1e_igkin_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_NUC)
            have = int1e_ignuc_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_RXP)
            have = int1e_giao_irjxp_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_A11)
            have = int1e_giao_a11part_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_A01GP)
            have = int1e_a01gp_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_PRINVXP)
            have = int1e_prinvxp_cart(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case default
            have = .false.
         end select
      else
         select case (which)
         case (GIAO_OVLP)
            have = int1e_igovlp_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_KIN)
            have = int1e_igkin_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_NUC)
            have = int1e_ignuc_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_RXP)
            have = int1e_giao_irjxp_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_A11)
            have = int1e_giao_a11part_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_A01GP)
            have = int1e_a01gp_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case (GIAO_PRINVXP)
            have = int1e_prinvxp_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         case default
            have = .false.
         end select
      end if
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
      ! and only the symmetry of the second pair used, so this is the full
      ! fourth power of the basis size and the slowest step of a shielding.
      ! A bound on the differentiated integral would let it screen as the
      ! Fock build does.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! Symmetric
      real(dp), allocatable, intent(out) :: vj(:, :, :), vk(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: error

      integer, allocatable, target :: atm(:), bas(:)
      real(dp), allocatable, target :: env(:)
      real(dp), allocatable :: buf(:), vj_t(:, :, :), vk_t(:, :, :)
      type(cint_ws) :: ws
      integer :: dims(0:3), shls(0:3)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo
      integer :: i, j, k, l, comp, mx, idx, n
      real(dp) :: val
      logical :: have

      if (error%has_error()) return
      n = mol%nao
      allocate (atm(size(mol%atm)), bas(size(mol%bas)))
      atm = reshape(mol%atm, [size(mol%atm)])
      bas = reshape(mol%bas, [size(mol%bas)])
      allocate (env(size(mol%env)), source=mol%env)
      mx = 0
      do ish = 1, mol%nbas
         mx = max(mx, shell_dim(mol%cartesian, ish - 1, mol%bas))
      end do
      allocate (vj(n, n, 3), vk(n, n, 3))
      vj = 0.0_dp
      vk = 0.0_dp

      !$omp parallel default(none) shared(mol, density, vj, vk, atm, bas, env, mx, n) &
      !$omp    private(ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, i, j, k, l, &
      !$omp            comp, idx, val, shls, dims, buf, ws, have, vj_t, vk_t)
      ws%opt => null()
      allocate (buf(3*mx**4), vj_t(n, n, 3), vk_t(n, n, 3))
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
                  dims = [di, dj, dk, dl]
                  if (mol%cartesian) then
                     have = int2e_ig1_cart(buf(1:3*di*dj*dk*dl), dims, shls, atm, &
                                           mol%natm, bas, mol%nbas, env, ws)
                  else
                     have = int2e_ig1_sph(buf(1:3*di*dj*dk*dl), dims, shls, atm, &
                                          mol%natm, bas, mol%nbas, env, ws)
                  end if
                  if (.not. have) cycle
                  idx = 0
                  do comp = 1, 3
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
