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
   use libcint_fortran, only: LIBCINT_PTR_RINV_ORIG, LIBCINT_ATM_SLOTS, LIBCINT_ATOM_OF, &
                              LIBCINT_CHARGE_OF, LIBCINT_PTR_COORD, LIBCINT_NUC_MOD_OF, &
                              LIBCINT_PTR_ZETA, &
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
   public :: giao_smeared_charges
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
   real(dp), parameter :: GIAO_SCREEN_TOL = 1.0e-12_dp
      !! Estimated contribution to a first-order Fock element below which a
      !! quartet of the two-electron field derivative is skipped
   real(dp), parameter :: FUNCTION_EXTENT = 5.0_dp
      !! Bohr added to a centre's distance from the origin when bounding how
      !! large the gauge factor is over a pair of functions

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

   subroutine giao_smeared_charges(mol, coords, zeta, charges, matrix, error)
      !! The field derivative of the potential of a set of Gaussian charges
      !!
      !! For the one-electron operator
      !! `-sum_k q_k erf(zeta_k |r - r_k|) / |r - r_k|`, which is what a
      !! continuum's surface charges add to the Fock matrix, the matrix of
      !! its derivative through the gauge factor.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)   !! (3, n_charges), Bohr
      real(dp), intent(in) :: zeta(:)        !! The `zeta_k` of `erf(zeta_k r)`
      real(dp), intent(in) :: charges(:)
      real(dp), allocatable, intent(out) :: matrix(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: error

      ! libcint's code for a nucleus that is a Gaussian rather than a point,
      ! which `libcint_fortran` has no name for.
      integer, parameter :: GAUSSIAN_NUCLEUS = 2
      ! Slots added to `env` for the charge: three coordinates, the width,
      ! and the fractional charge libfint counts one past the width.
      integer, parameter :: N_CHARGE_SLOTS = 5
      integer, allocatable :: atm(:, :), bas(:, :)
      real(dp), allocatable :: env(:), buf(:), acc(:)
      integer :: shls(2)
      integer :: n_charges, ish, jsh, di, dj, io, jo, i, j, comp, mx, k, slot, n_env
      logical :: have

      if (error%has_error()) return
      n_charges = size(charges)
      if (size(coords, 2) /= n_charges .or. size(zeta) /= n_charges) then
         call error%set(ERROR_VALIDATION, "the smeared charges, their positions and "// &
                        "their widths differ in number")
         return
      end if
      allocate (matrix(mol%nao, mol%nao, N_VECTOR))
      matrix = 0.0_dp
      if (n_charges == 0) return

      ! There is no integral for one such charge at an arbitrary point, and
      ! the nuclear one sums over every atom it is told of. So the charge is
      ! put first in the atom list as a Gaussian nucleus of unit charge, the
      ! molecule's atoms after it with every shell's atom index moved up by
      ! one, and the integral told there is one atom: it sums over the charge
      ! alone and still finds every shell's centre. That reads `atm` past the
      ! one atom declared, which holds because the array is passed whole.
      n_env = size(mol%env)
      allocate (atm(LIBCINT_ATM_SLOTS, mol%natm + 1))
      allocate (bas, source=mol%bas)
      atm(:, 1) = 0
      atm(:, 2:) = mol%atm
      atm(LIBCINT_CHARGE_OF, 1) = 1
      atm(LIBCINT_PTR_COORD, 1) = n_env
      atm(LIBCINT_NUC_MOD_OF, 1) = GAUSSIAN_NUCLEUS
      atm(LIBCINT_PTR_ZETA, 1) = n_env + 3
      bas(LIBCINT_ATOM_OF, :) = bas(LIBCINT_ATOM_OF, :) + 1
      mx = max_shell_dim(mol)
      slot = n_env + 1

      ! The gauge factor vanishes between two functions on one centre and
      ! changes sign when they are exchanged, so only pairs of shells on
      ! different atoms are computed, each once.
      !$omp parallel default(none) &
      !$omp    shared(mol, matrix, atm, bas, mx, coords, zeta, charges, n_charges, n_env, slot) &
      !$omp    private(ish, jsh, di, dj, io, jo, i, j, comp, shls, buf, acc, have, env, k)
      allocate (buf(mx*mx*N_VECTOR), acc(mx*mx*N_VECTOR), env(n_env + N_CHARGE_SLOTS))
      env(1:n_env) = mol%env
      env(n_env + 1:) = 0.0_dp
      !$omp do schedule(dynamic)
      do ish = 2, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, ish - 1
            if (mol%bas(LIBCINT_ATOM_OF, ish) == mol%bas(LIBCINT_ATOM_OF, jsh)) cycle
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]
            acc(1:di*dj*N_VECTOR) = 0.0_dp
            do k = 1, n_charges
               env(slot:slot + 2) = coords(:, k)
               ! A Gaussian nucleus of exponent `a` has the potential
               ! `erf(sqrt(a) r) / r`.
               env(slot + 3) = zeta(k)**2
               if (mol%cartesian) then
                  have = libcint_1e_ignuc_cart(buf, shls, atm, 1, bas, mol%nbas, env) /= 0
               else
                  have = libcint_1e_ignuc_sph(buf, shls, atm, 1, bas, mol%nbas, env) /= 0
               end if
               if (have) acc(1:di*dj*N_VECTOR) = acc(1:di*dj*N_VECTOR) &
                                                 + charges(k)*buf(1:di*dj*N_VECTOR)
            end do
            do comp = 1, N_VECTOR
               do j = 1, dj
                  do i = 1, di
                     matrix(io + i, jo + j, comp) = acc(i + di*(j - 1 + dj*(comp - 1)))
                     matrix(jo + j, io + i, comp) = -acc(i + di*(j - 1 + dj*(comp - 1)))
                  end do
               end do
            end do
         end do
      end do
      !$omp end do
      deallocate (buf, acc, env)
      !$omp end parallel
   end subroutine giao_smeared_charges

   subroutine giao_coulomb_exchange(mol, density, vj, vk, error, bounds, coulomb_only)
      !! The field derivative of the two-electron integrals against a density
      !!
      !! With `(g i j|k l)` the integral differentiated through the gauge
      !! factor of the first pair,
      !!
      !!     vj_ij = sum_kl (g i j|k l) D_lk,      vk_il = sum_jk (g i j|k l) D_jk
      !!
      !! for each of the three field directions. Neither is symmetrised or
      !! signed here; the caller applies the convention it wants.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! Symmetric
      real(dp), allocatable, intent(out) :: vj(:, :, :), vk(:, :, :)   !! (n_ao, n_ao, 3)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: bounds(:, :)
         !! `schwarz_bounds` of `mol`, (n_shells, n_shells). With them a
         !! quartet whose estimated contribution is below `GIAO_SCREEN_TOL` is
         !! skipped; without them every quartet is computed.
      logical, intent(in), optional :: coulomb_only
         !! Leave `vk` zero: a reference with no exact exchange has no use for
         !! it, and more quartets can then be skipped. Off when absent.

      ! The gauge factor of a pair is `(1/2) (R_i - R_j) x r`, so the
      ! differentiated integral
      !   * vanishes when the two functions share a centre,
      !   * changes sign when they are exchanged, and
      !   * keeps the ordinary integral's symmetry in the second pair.
      ! Only pairs `i > j` on different atoms and `k >= l` are computed, a
      ! quarter of the quartets at most, and each lands in the four places
      ! the full sum would have put it.
      real(dp), allocatable :: buf(:), vj_t(:, :, :), vk_t(:, :, :)
      real(dp), allocatable :: dmax(:, :), reach(:, :)
      integer, allocatable :: pair_i(:), pair_j(:), shell_atom(:)
      integer :: shls(4)
      integer :: ret
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo
      integer :: i, j, k, l, comp, mx, idx, n, n_pairs, ipair, nsh
      real(dp) :: val, weight, estimate, d_here
      logical :: screen, with_exchange

      if (error%has_error()) return
      n = mol%nao
      nsh = mol%nbas
      with_exchange = .true.
      if (present(coulomb_only)) with_exchange = .not. coulomb_only
      mx = max_shell_dim(mol)
      allocate (vj(n, n, N_VECTOR), vk(n, n, N_VECTOR))
      vj = 0.0_dp
      vk = 0.0_dp
      screen = present(bounds)

      allocate (shell_atom(nsh))
      do ish = 1, nsh
         shell_atom(ish) = mol%bas(LIBCINT_ATOM_OF, ish) + 1
      end do

      ! The pairs that can be non-zero, and for each how far the gauge factor
      ! can reach: `(1/2) |R_i - R_j|` times the distance from the origin of
      ! the farther centre, padded by the functions' own extent.
      allocate (pair_i(nsh*(nsh - 1)/2), pair_j(nsh*(nsh - 1)/2), reach(nsh, nsh))
      reach = 0.0_dp
      n_pairs = 0
      do ish = 2, nsh
         do jsh = 1, ish - 1
            if (shell_atom(ish) == shell_atom(jsh)) cycle
            n_pairs = n_pairs + 1
            pair_i(n_pairs) = ish
            pair_j(n_pairs) = jsh
            reach(ish, jsh) = 0.5_dp*norm2(mol%coords(:, shell_atom(ish)) &
                                           - mol%coords(:, shell_atom(jsh))) &
                              *(max(norm2(mol%coords(:, shell_atom(ish))), &
                                    norm2(mol%coords(:, shell_atom(jsh)))) + FUNCTION_EXTENT)
         end do
      end do

      ! Largest density element between each pair of shells.
      allocate (dmax(nsh, nsh))
      do jsh = 1, nsh
         do ish = 1, nsh
            dmax(ish, jsh) = maxval(abs(density(mol%shell_offset(ish) + 1:mol%shell_offset(ish + 1), &
                                                mol%shell_offset(jsh) + 1:mol%shell_offset(jsh + 1))))
         end do
      end do

      !$omp parallel default(none) &
      !$omp    shared(mol, density, vj, vk, mx, n, nsh, n_pairs, pair_i, pair_j, reach, dmax, &
      !$omp           bounds, screen, with_exchange) &
      !$omp    private(ipair, ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, i, j, k, l, &
      !$omp            comp, idx, val, weight, estimate, d_here, shls, buf, ret, vj_t, vk_t)
      allocate (buf(N_VECTOR*mx**4), vj_t(n, n, N_VECTOR), vk_t(n, n, N_VECTOR))
      vj_t = 0.0_dp
      vk_t = 0.0_dp
      !$omp do schedule(dynamic)
      do ipair = 1, n_pairs
         ish = pair_i(ipair)
         jsh = pair_j(ipair)
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
         io = mol%shell_offset(ish)
         jo = mol%shell_offset(jsh)
         do ksh = 1, nsh
            dk = shell_dim(mol%cartesian, ksh - 1, mol%bas)
            ko = mol%shell_offset(ksh)
            do lsh = 1, ksh
               if (screen) then
                  ! The densities this quartet is contracted with: the second
                  ! pair's for the Coulomb term, the four cross blocks for
                  ! exchange.
                  d_here = dmax(ksh, lsh)
                  if (with_exchange) then
                     d_here = max(d_here, dmax(jsh, ksh), dmax(jsh, lsh), &
                                  dmax(ish, ksh), dmax(ish, lsh))
                  end if
                  estimate = reach(ish, jsh)*bounds(ish, jsh)*bounds(ksh, lsh)*d_here
                  if (estimate < GIAO_SCREEN_TOL) cycle
               end if
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
               ! A second pair of two different shells stands for both of
               ! its orders in the Coulomb sum.
               weight = 2.0_dp
               if (ksh == lsh) weight = 1.0_dp
               idx = 0
               if (.not. with_exchange) then
                  do comp = 1, N_VECTOR
                     do l = lo + 1, lo + dl
                        do k = ko + 1, ko + dk
                           do j = jo + 1, jo + dj
                              do i = io + 1, io + di
                                 idx = idx + 1
                                 vj_t(i, j, comp) = vj_t(i, j, comp) &
                                                    + weight*buf(idx)*density(l, k)
                              end do
                           end do
                        end do
                     end do
                  end do
                  cycle
               end if
               do comp = 1, N_VECTOR
                  do l = lo + 1, lo + dl
                     do k = ko + 1, ko + dk
                        do j = jo + 1, jo + dj
                           do i = io + 1, io + di
                              idx = idx + 1
                              val = buf(idx)
                              vj_t(i, j, comp) = vj_t(i, j, comp) + weight*val*density(l, k)
                              vk_t(i, l, comp) = vk_t(i, l, comp) + val*density(j, k)
                              vk_t(j, l, comp) = vk_t(j, l, comp) - val*density(i, k)
                              if (ksh /= lsh) then
                                 vk_t(i, k, comp) = vk_t(i, k, comp) + val*density(j, l)
                                 vk_t(j, k, comp) = vk_t(j, k, comp) - val*density(i, l)
                              end if
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

      ! The Coulomb matrix is antisymmetric in its own two indices, and only
      ! the half with `i > j` was accumulated.
      do comp = 1, N_VECTOR
         vj(:, :, comp) = vj(:, :, comp) - transpose(vj(:, :, comp))
      end do
   end subroutine giao_coulomb_exchange

end module mqc_czt_giao
