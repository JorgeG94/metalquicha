!! Analytic nuclear gradients on the CPU backend
module mqc_libcint_gradient
   !! The derivative of a converged SCF energy with respect to the nuclear
   !! coordinates, for a restricted or unrestricted reference.
   !!
   !! **Why a variational wavefunction makes this cheap.** The energy is
   !! stationary with respect to the orbitals, so their derivatives do not
   !! appear: nothing here solves a response equation. What is left is the
   !! derivative of the *integrals*, contracted with densities the SCF already
   !! produced. That is why Hartree-Fock gradients come before MP2 or coupled
   !! cluster ones, which are not stationary and do need a response solve.
   !!
   !! **The four terms.**
   !!
   !!   * nuclear repulsion, which is analytic and needs no integrals
   !!   * the core Hamiltonian, `Tr(D dH/dR)`
   !!   * the two-electron term, `Tr(D dG/dR)`
   !!   * the Pulay term, `-Tr(W dS/dR)`, from the basis functions moving with
   !!     their atoms. This is the one that has no analogue in a plane-wave
   !!     code, and `W` is the energy-weighted density rather than the density
   !!
   !! **libcint differentiates the bra, and only the bra.** Every `ip` entry
   !! point returns the derivative with respect to the first shell's centre;
   !! the ket's derivative is recovered from translational invariance rather
   !! than integrated for. That is what the transposes and the factors of two
   !! below are doing, and it is why the sum of the returned gradient over
   !! atoms is zero by construction only if they are right -- which is worth
   !! checking, and is checked in the tests.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim, atom_ao_blocks
   use libcint_fortran, only: libcint_1e_ipovlp_sph, libcint_1e_ipovlp_cart, &
                              libcint_1e_ipkin_sph, libcint_1e_ipkin_cart, &
                              libcint_1e_ipnuc_sph, libcint_1e_ipnuc_cart, &
                              libcint_1e_iprinv_sph, libcint_1e_iprinv_cart, &
                              libcint_2e_ip1_sph, libcint_2e_ip1_cart, &
                              libcint_2e_ip1_sph_optimizer, libcint_2e_ip1_cart_optimizer, &
                              libcint_del_optimizer, LIBCINT_PTR_RINV_ORIG
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr
   implicit none
   private

   public :: libcint_scf_gradient
   public :: nuclear_repulsion_gradient   !! Exposed for the tests

   ! Which one-electron derivative `one_electron_deriv` should build.
   integer, parameter :: DERIV_OVLP = 1
   integer, parameter :: DERIV_KIN = 2
   integer, parameter :: DERIV_NUC = 3

contains

   subroutine libcint_scf_gradient(mol, density, density_beta, orbitals, orbitals_beta, &
                                   orbital_energies, orbital_energies_beta, &
                                   n_occupied, n_occupied_beta, gradient, error)
      !! dE/dR for a converged SCF, in Hartree/Bohr
      !!
      !! The beta arguments are what makes this one routine rather than two.
      !! Present means unrestricted: the Coulomb field is built from the total
      !! density and exchange from each spin separately, which is the only
      !! place the two paths differ. Absent means closed shell, where the
      !! density already carries its factor of two.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)          !! Alpha, or the total for RHF
      real(dp), intent(in), optional :: density_beta(:, :)
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in), optional :: orbitals_beta(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      real(dp), intent(in), optional :: orbital_energies_beta(:)
      integer, intent(in) :: n_occupied
      integer, intent(in), optional :: n_occupied_beta
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, natm)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: total_density(:, :), weighted(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :)
      real(dp), allocatable :: vhf(:, :, :), vhf_beta(:, :, :)
      real(dp), allocatable :: vrinv(:, :, :), hcore_a(:, :, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: nao, iatom, comp, p0, p1
      logical :: unrestricted

      unrestricted = present(density_beta)
      nao = mol%nao

      if (size(density, 1) /= nao) then
         call error%set(ERROR_VALIDATION, &
                        "libcint_scf_gradient: the density does not match this basis")
         return
      end if

      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp

      call nuclear_repulsion_gradient(mol, gradient)

      ! The total density drives the Coulomb field and the core Hamiltonian
      ! trace; exchange is per spin and is handled inside the two-electron
      ! contraction.
      allocate (total_density(nao, nao), source=density)
      if (unrestricted) total_density = total_density + density_beta

      call energy_weighted_density(orbitals, orbital_energies, n_occupied, &
                                   .not. unrestricted, weighted)
      if (unrestricted) then
         block
            real(dp), allocatable :: weighted_beta(:, :)
            call energy_weighted_density(orbitals_beta, orbital_energies_beta, &
                                         n_occupied_beta, .false., weighted_beta)
            weighted = weighted + weighted_beta
         end block
      end if

      ! The sign convention is libcint's: its `ip` integrals carry a nabla on
      ! the bra, and the derivative of the integral with respect to the atom
      ! the bra sits on is minus that.
      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)

      if (unrestricted) then
         call two_electron_deriv(mol, total_density, vhf, error, &
                                 exchange_density=density)
         if (error%has_error()) return
         call two_electron_deriv(mol, total_density, vhf_beta, error, &
                                 exchange_density=density_beta)
         if (error%has_error()) return
      else
         call two_electron_deriv(mol, total_density, vhf, error)
         if (error%has_error()) return
      end if

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      allocate (vrinv(nao, nao, 3), hcore_a(nao, nao, 3))

      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         ! dH/dR_A has two parts that look alike and are not. Moving atom A
         ! moves the basis functions centred on it, which is the h1 block
         ! below; it also moves the *nucleus*, and every electron feels that
         ! wherever its orbital sits. The second is the Hellmann-Feynman term,
         ! it involves no basis-set derivative at all, and it is what `iprinv`
         ! is for -- differentiating the origin of 1/|r-R| rather than a
         ! basis-function centre.
         call iprinv_deriv_at(mol, iatom, vrinv)
         vrinv = -mol%charges(iatom)*vrinv
         if (counts(iatom) > 0) then
            vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         end if
         do comp = 1, 3
            hcore_a(:, :, comp) = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
         end do

         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(hcore_a(:, :, comp)*total_density)
         end do

         if (counts(iatom) == 0) cycle

         ! Twice, because the nabla was on the bra: the ket's contribution is
         ! the same number by the symmetry of the integrals, so it is counted
         ! rather than computed.
         do comp = 1, 3
            if (unrestricted) then
               gradient(comp, iatom) = gradient(comp, iatom) &
                                       + 2.0_dp*sum(vhf(p0:p1, :, comp)*density(p0:p1, :)) &
                                       + 2.0_dp*sum(vhf_beta(p0:p1, :, comp)*density_beta(p0:p1, :))
            else
               gradient(comp, iatom) = gradient(comp, iatom) &
                                       + 2.0_dp*sum(vhf(p0:p1, :, comp)*total_density(p0:p1, :))
            end if
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*weighted(p0:p1, :))
         end do
      end do

   end subroutine libcint_scf_gradient

   subroutine nuclear_repulsion_gradient(mol, gradient)
      !! dE_NN/dR, accumulated into `gradient`
      !!
      !! E_NN = sum_{A<B} Z_A Z_B / R_AB, so the derivative on A points away
      !! from every other nucleus -- repulsion pushes atoms apart, and the
      !! gradient being the direction of steepest *ascent* makes the sign the
      !! one below.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(inout) :: gradient(:, :)

      integer :: ia, ib
      real(dp) :: rvec(3)
      real(dp) :: r

      do ia = 1, mol%natm
         do ib = 1, mol%natm
            if (ia == ib) cycle
            rvec = mol%coords(:, ia) - mol%coords(:, ib)
            r = norm2(rvec)
            gradient(:, ia) = gradient(:, ia) &
                              - mol%charges(ia)*mol%charges(ib)*rvec/r**3
         end do
      end do
   end subroutine nuclear_repulsion_gradient

   subroutine energy_weighted_density(orbitals, energies, n_occ, closed_shell, weighted)
      !! W = sum_i occ_i eps_i C_i C_i^T
      !!
      !! The Pulay term contracts against this rather than the density,
      !! because what moves with the basis functions is the orthonormality
      !! constraint, and its Lagrange multipliers are the orbital energies.
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: energies(:)
      integer, intent(in) :: n_occ
      logical, intent(in) :: closed_shell   !! Occupation two rather than one
      real(dp), allocatable, intent(out) :: weighted(:, :)

      integer :: nao, mu, nu, i
      real(dp) :: occupation, acc

      nao = size(orbitals, 1)
      allocate (weighted(nao, nao))
      weighted = 0.0_dp

      occupation = 1.0_dp
      if (closed_shell) occupation = 2.0_dp

      do nu = 1, nao
         do mu = 1, nao
            acc = 0.0_dp
            do i = 1, n_occ
               acc = acc + occupation*energies(i)*orbitals(mu, i)*orbitals(nu, i)
            end do
            weighted(mu, nu) = acc
         end do
      end do
   end subroutine energy_weighted_density

   subroutine one_electron_deriv(mol, matrix, which)
      !! A one-electron derivative matrix, as (nao, nao, 3)
      !!
      !! The same shell-pair loop as `one_electron` next door, with three
      !! components per block instead of one. libcint returns them with the
      !! component slowest, so a block is buf(di, dj, 1:3).
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), allocatable, intent(out) :: matrix(:, :, :)
      integer, intent(in) :: which

      real(dp), allocatable :: buf(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, comp, ret, io, jo, mx

      allocate (matrix(mol%nao, mol%nao, 3))
      matrix = 0.0_dp

      ! Flat, and indexed by hand below. libcint packs a block with the
      ! *actual* shell dimensions, so a rank-3 buffer declared at the largest
      ! shell has the wrong strides for every smaller one -- which is invisible
      ! in a basis of s functions and wrong everywhere else.
      mx = largest_shell(mol)
      allocate (buf(mx*mx*3))

      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]

            if (mol%cartesian) then
               select case (which)
               case (DERIV_OVLP)
                  ret = libcint_1e_ipovlp_cart(buf, shls, mol%atm, mol%natm, &
                                               mol%bas, mol%nbas, mol%env)
               case (DERIV_KIN)
                  ret = libcint_1e_ipkin_cart(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, mol%env)
               case default
                  ret = libcint_1e_ipnuc_cart(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, mol%env)
               end select
            else
               select case (which)
               case (DERIV_OVLP)
                  ret = libcint_1e_ipovlp_sph(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, mol%nbas, mol%env)
               case (DERIV_KIN)
                  ret = libcint_1e_ipkin_sph(buf, shls, mol%atm, mol%natm, &
                                             mol%bas, mol%nbas, mol%env)
               case default
                  ret = libcint_1e_ipnuc_sph(buf, shls, mol%atm, mol%natm, &
                                             mol%bas, mol%nbas, mol%env)
               end select
            end if

            if (ret == 0) cycle
            do comp = 1, 3
               do j = 1, dj
                  do i = 1, di
                     matrix(io + i, jo + j, comp) = buf(i + di*(j - 1 + dj*(comp - 1)))
                  end do
               end do
            end do
         end do
      end do
   end subroutine one_electron_deriv

   subroutine iprinv_deriv_at(mol, iatom, matrix)
      !! The 1/|r-R| derivative with the operator origin on atom `iatom`
      !!
      !! `env` is copied rather than modified, because the origin lives in it
      !! and the molecule is shared. One copy per atom of an array that is a
      !! few kilobytes is not worth avoiding; mutating the caller's molecule
      !! and putting it back would be.
      type(libcint_molecule_t), intent(in) :: mol
      integer, intent(in) :: iatom
      real(dp), intent(out) :: matrix(:, :, :)

      real(dp), allocatable :: env(:), buf(:)
      integer :: shls(2)
      integer :: ish, jsh, di, dj, i, j, comp, ret, io, jo, mx

      allocate (env(size(mol%env)), source=mol%env)
      env(LIBCINT_PTR_RINV_ORIG + 1:LIBCINT_PTR_RINV_ORIG + 3) = mol%coords(:, iatom)

      matrix = 0.0_dp
      mx = largest_shell(mol)
      allocate (buf(mx*mx*3))

      do ish = 1, mol%nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, mol%nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            shls = [ish - 1, jsh - 1]

            if (mol%cartesian) then
               ret = libcint_1e_iprinv_cart(buf, shls, mol%atm, mol%natm, &
                                            mol%bas, mol%nbas, env)
            else
               ret = libcint_1e_iprinv_sph(buf, shls, mol%atm, mol%natm, &
                                           mol%bas, mol%nbas, env)
            end if

            if (ret == 0) cycle
            do comp = 1, 3
               do j = 1, dj
                  do i = 1, di
                     matrix(io + i, jo + j, comp) = buf(i + di*(j - 1 + dj*(comp - 1)))
                  end do
               end do
            end do
         end do
      end do
   end subroutine iprinv_deriv_at

   subroutine two_electron_deriv(mol, density, vhf, error, exchange_density)
      !! `J - K/2` built from the differentiated ERIs, as (nao, nao, 3)
      !!
      !! `density` drives the Coulomb field and is the total density in both
      !! the restricted and unrestricted cases. `exchange_density` is the one
      !! exchange is built from: absent for a closed shell, where the total
      !! density carries its factor of two and `J - K/2` is the right
      !! combination; present for one spin of an unrestricted reference, where
      !! exchange takes that spin alone and the coefficient is one rather than
      !! a half.
      !!
      !! **One permutation survives, and it is used.** A differentiated quartet
      !! has none of the eightfold symmetry of an ordinary one: the nabla
      !! distinguishes the first index from the second and the bra from the
      !! ket. What is left untouched is the ket pair, `(∇i j|k l) = (∇i j|l k)`,
      !! so this walks `l <= k` and counts the off-diagonal pair twice -- once
      !! as itself and once transposed, which is not the same contraction and
      !! has to be written out rather than doubled.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), allocatable, intent(out) :: vhf(:, :, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: exchange_density(:, :)

      real(dp), allocatable :: buf(:), vj(:, :, :), vk(:, :, :)
      real(dp), allocatable :: vj_local(:, :, :), vk_local(:, :, :)
      real(dp), allocatable :: exchange_from(:, :)
      real(dp) :: g, k_scale
      type(c_ptr) :: opt
      integer :: shls(4)
      integer :: ish, jsh, ksh, lsh, di, dj, dk, dl
      integer :: io, jo, ko, lo, i, j, k, l, comp, ret, mx, idx, nao, nbas

      mx = largest_shell(mol)
      nao = mol%nao
      nbas = mol%nbas
      allocate (vj(nao, nao, 3), vk(nao, nao, 3))
      vj = 0.0_dp
      vk = 0.0_dp

      allocate (exchange_from(nao, nao))
      if (present(exchange_density)) then
         exchange_from = exchange_density
         k_scale = 1.0_dp
      else
         exchange_from = density
         k_scale = 0.5_dp
      end if

      opt = c_null_ptr
      if (mol%cartesian) then
         call libcint_2e_ip1_cart_optimizer(opt, mol%atm, mol%natm, mol%bas, mol%nbas, mol%env)
      else
         call libcint_2e_ip1_sph_optimizer(opt, mol%atm, mol%natm, mol%bas, mol%nbas, mol%env)
      end if

      ! Thread-local accumulators merged once, rather than atomics on a shared
      ! matrix: the inner update is two scattered writes per integral, and
      ! making those atomic costs more than the integral. The price is
      ! 2 * nao^2 * 3 doubles per thread, which is worth watching on a few
      ! thousand functions -- the same caveat the direct Fock build carries.
      !
      ! schedule(dynamic) because the l <= k triangle makes the work per ish
      ! uneven, and a static split leaves threads idle on the tail.
      !$omp parallel default(none) &
      !$omp    shared(mol, density, exchange_from, opt, vj, vk, mx, nao, nbas) &
      !$omp    private(ish, jsh, ksh, lsh, di, dj, dk, dl, io, jo, ko, lo, &
      !$omp            i, j, k, l, comp, ret, idx, g, shls, buf, vj_local, vk_local)
      allocate (buf(mx**4*3))
      allocate (vj_local(nao, nao, 3), vk_local(nao, nao, 3))
      vj_local = 0.0_dp
      vk_local = 0.0_dp

      !$omp do schedule(dynamic)
      do ish = 1, nbas
         di = shell_dim(mol%cartesian, ish - 1, mol%bas)
         io = mol%shell_offset(ish)
         do jsh = 1, nbas
            dj = shell_dim(mol%cartesian, jsh - 1, mol%bas)
            jo = mol%shell_offset(jsh)
            do ksh = 1, nbas
               dk = shell_dim(mol%cartesian, ksh - 1, mol%bas)
               ko = mol%shell_offset(ksh)
               do lsh = 1, ksh
                  dl = shell_dim(mol%cartesian, lsh - 1, mol%bas)
                  lo = mol%shell_offset(lsh)
                  shls = [ish - 1, jsh - 1, ksh - 1, lsh - 1]

                  if (mol%cartesian) then
                     ret = libcint_2e_ip1_cart(buf, shls, mol%atm, mol%natm, &
                                               mol%bas, nbas, mol%env, opt)
                  else
                     ret = libcint_2e_ip1_sph(buf, shls, mol%atm, mol%natm, &
                                              mol%bas, nbas, mol%env, opt)
                  end if
                  if (ret == 0) cycle

                  do comp = 1, 3
                     do l = 1, dl
                        do k = 1, dk
                           do j = 1, dj
                              do i = 1, di
                                 idx = i + di*(j - 1 + dj*(k - 1 + dk*(l - 1 + dl*(comp - 1))))
                                 g = buf(idx)
                                 vj_local(io + i, jo + j, comp) = vj_local(io + i, jo + j, comp) &
                                                                  + g*density(lo + l, ko + k)
                                 vk_local(io + i, lo + l, comp) = vk_local(io + i, lo + l, comp) &
                                                                  + g*exchange_from(jo + j, ko + k)
                                 ! The transposed ket pair. Skipped on the
                                 ! diagonal, where it is the same quartet.
                                 if (lsh /= ksh) then
                                    vj_local(io + i, jo + j, comp) = vj_local(io + i, jo + j, comp) &
                                                                     + g*density(ko + k, lo + l)
                                    vk_local(io + i, ko + k, comp) = vk_local(io + i, ko + k, comp) &
                                                                     + g*exchange_from(jo + j, lo + l)
                                 end if
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
      vj = vj + vj_local
      vk = vk + vk_local
      !$omp end critical
      deallocate (buf, vj_local, vk_local)
      !$omp end parallel

      call libcint_del_optimizer(opt)

      ! Minus, for the same reason the one-electron ones carry it: libcint
      ! returns the nabla on the bra, and the derivative with respect to the
      ! atom is its negative.
      allocate (vhf(nao, nao, 3))
      vhf = -(vj - k_scale*vk)

   end subroutine two_electron_deriv

   pure function largest_shell(mol) result(n)
      !! Functions in the biggest shell, for sizing a scratch block
      type(libcint_molecule_t), intent(in) :: mol
      integer :: n

      integer :: ish

      n = 1
      do ish = 1, mol%nbas
         n = max(n, mol%shell_offset(ish + 1) - mol%shell_offset(ish))
      end do
   end function largest_shell

end module mqc_libcint_gradient
