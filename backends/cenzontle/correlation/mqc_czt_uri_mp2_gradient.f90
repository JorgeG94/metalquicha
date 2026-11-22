!! The RI-MP2 nuclear gradient over an unrestricted reference
module mqc_czt_uri_mp2_gradient
   !! An exact unrestricted Hartree-Fock reference with a fitted correlation
   !! energy, the open-shell twin of `mqc_czt_ri_mp2_gradient`.
   !!
   !! Written from the Hylleraas functional, which is stationary in the
   !! amplitudes, so that only the orbitals respond:
   !!
   !!     E2 = sum_s [ 2 sum_iaP Gamma^P_ia B^P_ia  +  sum_pq P_pq f_pq ]
   !!          evaluated where it equals the MP2 energy
   !!
   !! with `B^P_ia` the fitted integrals of one spin, `f` that spin's Fock
   !! matrix in its own orbitals, `P` the unrelaxed correlation density and
   !! `Gamma` the three-index density
   !!
   !!     Gamma^P_ia = sum_jb T_ij^ab B^P_jb  +  sum_JB t_iJ^aB B^P_JB
   !!
   !! in which `T` is the antisymmetrised same-spin amplitude and `t` the
   !! opposite-spin one. Everything below is that functional differentiated:
   !! with respect to the integrals, directly, and with respect to the
   !! orbitals, through one coupled-perturbed solve.
   !!
   !! **Per spin, and no factor of two anywhere a closed shell has one.** The
   !! two spins meet in three places only: the opposite-spin amplitudes, the
   !! Coulomb part of every mean field, and the coupled-perturbed operator.
   !!
   !! **The reference is exact and so is the response.** The coupled-perturbed
   !! operator is `response_product_uhf`, integral-direct, which is the
   !! operator the reference's orbitals are stationary for. Fitting it, as the
   !! closed-shell gradient does by default, would differentiate an energy
   !! this is not.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_logger, only: logger => global_logger
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks, build_df_mo_tensor, &
                                three_centre, two_centre, metric_inverse_sqrt
   use mqc_czt_gradient, only: nuclear_repulsion_gradient, one_electron_deriv, &
                               nuclear_attraction_deriv_at, two_electron_gradient, &
                               DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_czt_direct, only: schwarz_bounds
   use mqc_czt_response_product, only: response_mean_field_uhf, response_product_uhf
   use mqc_czt_ri_mp2_gradient, only: fitted_two_particle_gradient
   implicit none
   private

   public :: czt_uri_mp2_gradient

   integer, parameter :: Z_MAX_ITER = 200
   real(dp), parameter :: Z_TOL = 1.0e-11_dp
      !! Largest element of the Z-vector equations' residual at convergence.
      !! The residual enters the gradient to first order, not second.
   real(dp), parameter :: SINGULAR_PIVOT = 1.0e-13_dp
      !! A pivot this small beside the projected operator's largest element
      !! means the operator is singular on the trial space.
   integer, parameter :: LINE_LEN = 160

   type :: spin_block_t
      !! One spin's orbitals, fitted integrals and densities
      real(dp), allocatable :: c(:, :)          !! (n_ao, n_mo)
      real(dp), allocatable :: eps(:)           !! (n_mo)
      integer :: n_o = 0                        !! Occupied, frozen included
      integer :: n_act = 0                      !! Occupied and correlated
      integer :: n_v = 0
      integer :: n_ov = 0                       !! `n_act*n_v`
      real(dp), allocatable :: b(:, :)
         !! (n_ov, n_aux) fitted integrals, the pair index `i + (a-1)*n_act`
      real(dp), allocatable :: half(:, :)
         !! (n_ov, n_aux) the three-index density before the metric, which is
         !! what the energy and the metric term contract with
      real(dp), allocatable :: gamma(:, :, :)
         !! (n_aux, n_act, n_v) the three-index density with the inverse
         !! square root of the metric applied
      real(dp), allocatable :: p(:, :)
         !! (n_mo, n_mo) the correlation density: unrelaxed occupied-occupied
         !! and virtual-virtual blocks first, relaxed blocks added later
      real(dp), allocatable :: l1(:, :)
         !! (n_mo, n_mo) twice the three-index density against `(pq|P)`: row
         !! any orbital, column a correlated occupied or a virtual one. Zero
         !! in frozen columns.
      real(dp), allocatable :: w(:, :)          !! (n_mo, n_mo) energy-weighted density
   end type spin_block_t

contains

   subroutine czt_uri_mp2_gradient(mol, aux, coeff_a, coeff_b, eps_a, eps_b, n_occ_a, &
                                   n_occ_b, gradient, error, n_frozen, e_corr, &
                                   relaxed_alpha, relaxed_beta)
      !! dE(UHF + RI-MP2)/dR for an unrestricted reference, in Hartree/Bohr
      !!
      !! The whole gradient, reference included: its one- and two-electron
      !! terms are contracted with the reference density plus the relaxed
      !! correlation density in one pass each.
      type(czt_molecule_t), intent(in) :: mol
      type(czt_molecule_t), intent(in) :: aux   !! Auxiliary basis, same atoms
      real(dp), intent(in) :: coeff_a(:, :), coeff_b(:, :)   !! (n_ao, n_mo), canonical
      real(dp), intent(in) :: eps_a(:), eps_b(:)
      integer, intent(in) :: n_occ_a, n_occ_b
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, n_atoms)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen
         !! Core orbitals left out of the correlation, the same count in each
         !! spin, as the energy takes it
      real(dp), intent(out), optional :: e_corr
         !! The correlation energy recomputed from the three-index density.
         !! It is the energy routine's number by another route, for a caller
         !! that wants to know the two agree.
      real(dp), allocatable, intent(out), optional :: relaxed_alpha(:, :), relaxed_beta(:, :)
         !! (n_ao, n_ao) each spin's relaxed density, reference included

      type(spin_block_t), target :: sp(2)
      real(dp), allocatable :: three_ao(:, :), metric(:, :), jm12(:, :), bounds(:, :)
      real(dp), allocatable :: zero_h(:, :)
      real(dp), allocatable :: dens_a(:, :, :), dens_b(:, :, :), g_a(:, :, :), g_b(:, :, :)
      real(dp), allocatable :: x_a(:, :), x_b(:, :), z_a(:, :), z_b(:, :)
      real(dp), allocatable :: gamma_ao(:, :, :), gamma_pq(:, :)
      real(dp), allocatable :: d_ref(:, :, :), d_rel(:, :, :), w_ao(:, :), cross(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :)
      real(dp), allocatable :: vrinv(:, :, :), hcore(:, :), d_all(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, n_mo, n_aux, frozen, s, i, a, f, iatom, comp, p0, p1
      real(dp) :: e2

      if (error%has_error()) return
      n_ao = mol%nao
      n_mo = size(coeff_a, 2)
      frozen = 0
      if (present(n_frozen)) frozen = n_frozen
      if (size(coeff_b, 2) /= n_mo .or. size(eps_a) /= n_mo .or. size(eps_b) /= n_mo) then
         call error%set(ERROR_VALIDATION, "the unrestricted RI-MP2 gradient was given "// &
                        "orbitals and energies of different sizes")
         return
      end if
      if (frozen < 0 .or. frozen >= n_occ_a .or. n_occ_a >= n_mo) then
         call error%set(ERROR_VALIDATION, "the unrestricted RI-MP2 gradient needs a "// &
                        "correlated occupied orbital and a virtual one of alpha spin")
         return
      end if
      if (frozen > n_occ_b) then
         call error%set(ERROR_VALIDATION, "the unrestricted RI-MP2 gradient's frozen "// &
                        "core is larger than the number of beta electrons")
         return
      end if
      if (aux%natm /= mol%natm) then
         call error%set(ERROR_VALIDATION, "the unrestricted RI-MP2 gradient's auxiliary "// &
                        "basis is on a different number of atoms than the orbital basis")
         return
      end if
      if (maxval(abs(aux%coords - mol%coords)) > 1.0e-10_dp) then
         call error%set(ERROR_VALIDATION, "the unrestricted RI-MP2 gradient's auxiliary "// &
                        "basis is on a different geometry than the orbital basis")
         return
      end if

      sp(1)%c = coeff_a
      sp(1)%eps = eps_a
      sp(1)%n_o = n_occ_a
      sp(2)%c = coeff_b
      sp(2)%eps = eps_b
      sp(2)%n_o = n_occ_b
      do s = 1, 2
         sp(s)%n_act = sp(s)%n_o - frozen
         sp(s)%n_v = n_mo - sp(s)%n_o
         sp(s)%n_ov = sp(s)%n_act*sp(s)%n_v
         allocate (sp(s)%p(n_mo, n_mo), sp(s)%l1(n_mo, n_mo), sp(s)%w(n_mo, n_mo))
         sp(s)%p = 0.0_dp
         sp(s)%l1 = 0.0_dp
         sp(s)%w = 0.0_dp
      end do

      ! ---- fitted integrals, amplitudes, and the two unrelaxed densities ---
      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, jm12, error)
      if (error%has_error()) return
      n_aux = size(jm12, 1)
      do s = 1, 2
         call fitted_block(mol, aux, sp(s), frozen, n_aux, error)
         if (error%has_error()) return
      end do
      call amplitudes_and_densities(sp, frozen, n_aux, e2)
      if (present(e_corr)) e_corr = e2
      do s = 1, 2
         call apply_metric(sp(s), jm12, n_aux)
      end do

      ! ---- the three-index density against (pq|P) --------------------------
      call three_centre(mol, aux, three_ao)
      do s = 1, 2
         call lagrangian_from_gamma(sp(s), three_ao, frozen, n_ao, n_mo, n_aux)
      end do
      deallocate (three_ao)

      ! ---- the occupied-frozen rotations -----------------------------------
      !
      ! The amplitudes never touch a frozen orbital, but the energy is not
      ! invariant to mixing one with a correlated one, so that rotation
      ! responds. Both orbitals are occupied, the reference density does not
      ! change, and the equation is diagonal: the Lagrangian's element over
      ! the orbital-energy gap. Stored negated and halved on each side, as
      ! every relaxed block of `p` is, so `p` contracts with a Fock
      ! derivative as it stands.
      do s = 1, 2
         do i = frozen + 1, sp(s)%n_o
            do f = 1, frozen
               sp(s)%p(f, i) = -0.5_dp*sp(s)%l1(f, i)/(sp(s)%eps(f) - sp(s)%eps(i))
               sp(s)%p(i, f) = sp(s)%p(f, i)
            end do
         end do
      end do

      ! ---- the Z-vector ----------------------------------------------------
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return
      allocate (zero_h(n_ao, n_ao), dens_a(n_ao, n_ao, 1), dens_b(n_ao, n_ao, 1))
      zero_h = 0.0_dp
      dens_a(:, :, 1) = to_ao(sp(1)%c, sp(1)%p)
      dens_b(:, :, 1) = to_ao(sp(2)%c, sp(2)%p)
      call response_mean_field_uhf(mol, dens_a, dens_b, zero_h, bounds, g_a, g_b, error)
      if (error%has_error()) return

      ! Right-hand side, virtual by occupied: the Lagrangian's antisymmetric
      ! part, and the reference's mean field of the density so far. A frozen
      ! column has no `l1(a, i)` of its own and still has the other two.
      call z_rhs(sp(1), g_a(:, :, 1), x_a)
      call z_rhs(sp(2), g_b(:, :, 1), x_b)
      call solve_z(mol, sp, zero_h, bounds, x_a, x_b, z_a, z_b, error)
      if (error%has_error()) return
      do i = 1, sp(1)%n_o
         do a = 1, sp(1)%n_v
            sp(1)%p(sp(1)%n_o + a, i) = -0.5_dp*z_a(a, i)
            sp(1)%p(i, sp(1)%n_o + a) = -0.5_dp*z_a(a, i)
         end do
      end do
      do i = 1, sp(2)%n_o
         do a = 1, sp(2)%n_v
            sp(2)%p(sp(2)%n_o + a, i) = -0.5_dp*z_b(a, i)
            sp(2)%p(i, sp(2)%n_o + a) = -0.5_dp*z_b(a, i)
         end do
      end do

      ! ---- the energy-weighted density -------------------------------------
      !
      ! The reference's mean field of the *relaxed* density, between occupied
      ! orbitals: the overlap derivative moves the reference density, and the
      ! Fock matrix the correlation density is contracted with moves with it.
      dens_a(:, :, 1) = to_ao(sp(1)%c, sp(1)%p)
      dens_b(:, :, 1) = to_ao(sp(2)%c, sp(2)%p)
      deallocate (g_a, g_b)
      call response_mean_field_uhf(mol, dens_a, dens_b, zero_h, bounds, g_a, g_b, error)
      if (error%has_error()) return
      call weighted_density(sp(1), g_a(:, :, 1), z_a, frozen)
      call weighted_density(sp(2), g_b(:, :, 1), z_b, frozen)

      ! ---- assembly --------------------------------------------------------
      allocate (d_ref(n_ao, n_ao, 2), d_rel(n_ao, n_ao, 2), w_ao(n_ao, n_ao))
      w_ao = 0.0_dp
      do s = 1, 2
         d_ref(:, :, s) = matmul(sp(s)%c(:, 1:sp(s)%n_o), &
                                 transpose(sp(s)%c(:, 1:sp(s)%n_o)))
         d_rel(:, :, s) = to_ao(sp(s)%c, sp(s)%p)
         w_ao = w_ao + to_ao(sp(s)%c, sp(s)%w)
      end do
      if (present(relaxed_alpha)) relaxed_alpha = d_ref(:, :, 1) + d_rel(:, :, 1)
      if (present(relaxed_beta)) relaxed_beta = d_ref(:, :, 2) + d_rel(:, :, 2)

      allocate (gradient(3, mol%natm), cross(3, mol%natm))
      gradient = 0.0_dp
      call nuclear_repulsion_gradient(mol, gradient)

      ! The reference's two-electron energy is quadratic in the densities, so
      ! its gradient at reference-plus-relaxed, less its gradient at relaxed
      ! alone, is the reference's own term and the cross term together.
      call two_electron_gradient(mol, sum(d_ref, 3) + sum(d_rel, 3), gradient, error, &
                                 density_alpha=d_ref(:, :, 1) + d_rel(:, :, 1), &
                                 density_beta=d_ref(:, :, 2) + d_rel(:, :, 2))
      if (error%has_error()) return
      cross = 0.0_dp
      call two_electron_gradient(mol, sum(d_rel, 3), cross, error, &
                                 density_alpha=d_rel(:, :, 1), density_beta=d_rel(:, :, 2))
      if (error%has_error()) return
      gradient = gradient - cross

      ! The fitted two-particle term. `fitted_two_particle_gradient` was
      ! written for a closed shell's density, which is twice a spin's, so
      ! both of its inputs are halved here.
      call fitted_densities(sp, jm12, n_ao, n_aux, gamma_ao, gamma_pq)
      call fitted_two_particle_gradient(mol, aux, gamma_ao, gamma_pq, gradient)

      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (vrinv(n_ao, n_ao, 3), hcore(n_ao, n_ao), d_all(n_ao, n_ao))
      d_all = sum(d_ref, 3) + sum(d_rel, 3)
      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)
         call nuclear_attraction_deriv_at(mol, iatom, vrinv)
         if (counts(iatom) > 0) vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         do comp = 1, 3
            hcore = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
            gradient(comp, iatom) = gradient(comp, iatom) + sum(hcore*d_all)
         end do
         if (counts(iatom) == 0) cycle
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*w_ao(p0:p1, :))
         end do
      end do
   end subroutine czt_uri_mp2_gradient

   function to_ao(c, mo) result(ao)
      !! `C M C^T`
      real(dp), intent(in) :: c(:, :), mo(:, :)
      real(dp) :: ao(size(c, 1), size(c, 1))
      real(dp) :: work(size(c, 1), size(c, 2))

      call pic_gemm(c, mo, work, beta=0.0_dp)
      call pic_gemm(work, c, ao, transb="T", beta=0.0_dp)
   end function to_ao

   subroutine fitted_block(mol, aux, spin, frozen, n_aux, error)
      !! `B^P_ia` for one spin, in the compound pair index
      type(czt_molecule_t), intent(in) :: mol, aux
      type(spin_block_t), intent(inout) :: spin
      integer, intent(in) :: frozen, n_aux
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: raw(:, :, :)
      integer :: i, p

      allocate (spin%b(spin%n_ov, n_aux))
      if (spin%n_ov == 0) return
      ! `build_df_mo_tensor` lays its result out `(a, P, i)`.
      call build_df_mo_tensor(mol, aux, spin%c(:, frozen + 1:spin%n_o), &
                              spin%c(:, spin%n_o + 1:), raw, error)
      if (error%has_error()) return
      do p = 1, n_aux
         do i = 1, spin%n_act
            spin%b(i:spin%n_ov:spin%n_act, p) = raw(:, p, i)
         end do
      end do
   end subroutine fitted_block

   subroutine amplitudes_and_densities(sp, frozen, n_aux, e2)
      !! First-order amplitudes, the three-index density, and the unrelaxed
      !! occupied-occupied and virtual-virtual correlation densities
      !!
      !! With `N = (1/2) T T^T + t t^T` in one spin's pair index, `T` the
      !! antisymmetrised same-spin amplitudes and `t` the opposite-spin ones:
      !!
      !!     P_ij = -sum_a N(ia, ja),      P_ab = +sum_i N(ia, ib)
      type(spin_block_t), intent(inout) :: sp(2)
      integer, intent(in) :: frozen, n_aux
      real(dp), intent(out) :: e2

      real(dp), allocatable :: same(:, :), opp(:, :), n(:, :)
      integer :: s, o

      e2 = 0.0_dp
      do s = 1, 2
         allocate (sp(s)%half(sp(s)%n_ov, n_aux))
         sp(s)%half = 0.0_dp
      end do
      ! Opposite spin first: one rectangle serves both spins.
      if (sp(1)%n_ov > 0 .and. sp(2)%n_ov > 0) then
         allocate (opp(sp(1)%n_ov, sp(2)%n_ov))
         call pic_gemm(sp(1)%b, sp(2)%b, opp, transb="T", beta=0.0_dp)
         call divide_by_gaps(opp, sp(1), sp(2), frozen)
         call pic_gemm(opp, sp(2)%b, sp(1)%half, beta=0.0_dp)
         call pic_gemm(opp, sp(1)%b, sp(2)%half, transa="T", beta=0.0_dp)
      end if
      do s = 1, 2
         if (sp(s)%n_ov == 0) cycle
         o = 3 - s
         allocate (same(sp(s)%n_ov, sp(s)%n_ov), n(sp(s)%n_ov, sp(s)%n_ov))
         call pic_gemm(sp(s)%b, sp(s)%b, same, transb="T", beta=0.0_dp)
         call antisymmetrize(same, sp(s)%n_act, sp(s)%n_v)
         call divide_by_gaps(same, sp(s), sp(s), frozen)
         call pic_gemm(same, sp(s)%b, sp(s)%half, beta=1.0_dp)
         call pic_gemm(same, same, n, transb="T", alpha=0.5_dp, beta=0.0_dp)
         if (allocated(opp)) then
            if (s == 1) then
               call pic_gemm(opp, opp, n, transb="T", beta=1.0_dp)
            else
               call pic_gemm(opp, opp, n, transa="T", beta=1.0_dp)
            end if
         end if
         call fold_pair_matrix(n, sp(s), frozen)
         deallocate (same, n)
         ! Each spin's half of `sum Gamma B`, which is twice the energy.
         e2 = e2 + 0.5_dp*sum(sp(s)%half*sp(s)%b)
      end do
   end subroutine amplitudes_and_densities

   subroutine antisymmetrize(g, n_act, n_v)
      !! `(ia|jb) - (ib|ja)` in place, in the pair index
      real(dp), intent(inout) :: g(:, :)
      integer, intent(in) :: n_act, n_v

      real(dp), allocatable :: copy(:, :)
      integer :: i, j, a, b

      copy = g
      do b = 1, n_v
         do j = 1, n_act
            do a = 1, n_v
               do i = 1, n_act
                  g(i + (a - 1)*n_act, j + (b - 1)*n_act) = &
                     copy(i + (a - 1)*n_act, j + (b - 1)*n_act) &
                     - copy(i + (b - 1)*n_act, j + (a - 1)*n_act)
               end do
            end do
         end do
      end do
   end subroutine antisymmetrize

   subroutine divide_by_gaps(g, left, right, frozen)
      !! Integrals to amplitudes: `g(ia, jb)/(e_i + e_j - e_a - e_b)`
      real(dp), intent(inout) :: g(:, :)
      type(spin_block_t), intent(in) :: left, right   !! The spins of `ia` and of `jb`
      integer, intent(in) :: frozen

      integer :: i, j, a, b

      do b = 1, right%n_v
         do j = 1, right%n_act
            do a = 1, left%n_v
               do i = 1, left%n_act
                  g(i + (a - 1)*left%n_act, j + (b - 1)*right%n_act) = &
                     g(i + (a - 1)*left%n_act, j + (b - 1)*right%n_act) &
                     /(left%eps(frozen + i) + right%eps(frozen + j) &
                       - left%eps(left%n_o + a) - right%eps(right%n_o + b))
               end do
            end do
         end do
      end do
   end subroutine divide_by_gaps

   subroutine fold_pair_matrix(n, spin, frozen)
      !! The pair-index matrix `N` traced down to the two one-particle blocks
      real(dp), intent(in) :: n(:, :)
      type(spin_block_t), intent(inout) :: spin
      integer, intent(in) :: frozen

      integer :: i, j, a, b, na

      na = spin%n_act
      do a = 1, spin%n_v
         do j = 1, na
            do i = 1, na
               spin%p(frozen + i, frozen + j) = spin%p(frozen + i, frozen + j) &
                                                - n(i + (a - 1)*na, j + (a - 1)*na)
            end do
         end do
      end do
      do b = 1, spin%n_v
         do a = 1, spin%n_v
            do i = 1, na
               spin%p(spin%n_o + a, spin%n_o + b) = spin%p(spin%n_o + a, spin%n_o + b) &
                                                    + n(i + (a - 1)*na, i + (b - 1)*na)
            end do
         end do
      end do
   end subroutine fold_pair_matrix

   subroutine apply_metric(spin, jm12, n_aux)
      !! The three-index density with `J^-1/2` applied, laid out `(P, i, a)`
      type(spin_block_t), intent(inout) :: spin
      real(dp), intent(in) :: jm12(:, :)
      integer, intent(in) :: n_aux

      real(dp), allocatable :: gm(:, :)
      integer :: i, a

      allocate (spin%gamma(n_aux, spin%n_act, spin%n_v))
      if (spin%n_ov == 0) return
      allocate (gm(spin%n_ov, n_aux))
      call pic_gemm(spin%half, jm12, gm, transb="T", beta=0.0_dp)
      do a = 1, spin%n_v
         do i = 1, spin%n_act
            spin%gamma(:, i, a) = gm(i + (a - 1)*spin%n_act, :)
         end do
      end do
   end subroutine apply_metric

   subroutine lagrangian_from_gamma(spin, three_ao, frozen, n_ao, n_mo, n_aux)
      !! `l1(p, i) = 2 sum_aP Gamma^P_ia (pa|P)` and
      !! `l1(p, a) = 2 sum_iP Gamma^P_ia (ip|P)`, `p` over every orbital
      type(spin_block_t), intent(inout) :: spin
      real(dp), intent(in) :: three_ao(:, :)   !! ((mu nu), P), raw
      integer, intent(in) :: frozen, n_ao, n_mo, n_aux

      real(dp), allocatable :: pq(:, :, :), block(:, :), slab(:, :), gslab(:, :)
      integer :: pp, i, a

      if (spin%n_ov == 0) return
      allocate (pq(n_mo, n_mo, n_aux), block(n_ao, n_ao))
      !$omp parallel do default(none) shared(pq, three_ao, spin, n_ao, n_aux) &
      !$omp    private(pp, block)
      do pp = 1, n_aux
         block = reshape(three_ao(:, pp), [n_ao, n_ao])
         pq(:, :, pp) = matmul(transpose(spin%c), matmul(block, spin%c))
      end do
      !$omp end parallel do

      allocate (slab(n_mo, n_aux))
      do a = 1, spin%n_v
         slab = pq(:, spin%n_o + a, :)
         call pic_gemm(slab, spin%gamma(:, :, a), &
                       spin%l1(:, frozen + 1:spin%n_o), alpha=2.0_dp, beta=1.0_dp)
      end do
      allocate (gslab(n_aux, spin%n_v))
      do i = 1, spin%n_act
         slab = pq(frozen + i, :, :)
         gslab = spin%gamma(:, i, :)
         call pic_gemm(slab, gslab, spin%l1(:, spin%n_o + 1:), alpha=2.0_dp, beta=1.0_dp)
      end do
   end subroutine lagrangian_from_gamma

   subroutine z_rhs(spin, g, x)
      !! `X_ai = l1(a, i) - l1(i, a) + 2 (C_v^T G C_o)_ai`
      type(spin_block_t), intent(in) :: spin
      real(dp), intent(in) :: g(:, :)   !! This spin's mean field of the density so far
      real(dp), allocatable, intent(out) :: x(:, :)   !! (n_v, n_o)

      real(dp), allocatable :: work(:, :)
      integer :: i, a, n_o

      n_o = spin%n_o
      allocate (x(spin%n_v, n_o))
      if (spin%n_v == 0 .or. n_o == 0) return
      allocate (work(size(g, 1), n_o))
      call pic_gemm(g, spin%c(:, 1:n_o), work, beta=0.0_dp)
      call pic_gemm(spin%c(:, n_o + 1:), work, x, transa="T", alpha=2.0_dp, beta=0.0_dp)
      do i = 1, n_o
         do a = 1, spin%n_v
            x(a, i) = x(a, i) + spin%l1(n_o + a, i) - spin%l1(i, n_o + a)
         end do
      end do
   end subroutine z_rhs

   subroutine solve_z(mol, sp, zero_h, bounds, x_a, x_b, z_a, z_b, error)
      !! `(A + B) z = X` over both spins, in a growing Krylov subspace
      !!
      !! The operator is the unrestricted orbital Hessian: symmetric always,
      !! and positive definite only where the reference is a minimum. A
      !! gradient is defined at any stationary point, so this is a Galerkin
      !! solve -- the equations projected onto the trial vectors and solved
      !! there -- which asks for a non-singular operator and not a definite
      !! one, where conjugate gradients would divide by a curvature that can
      !! be negative. Each new direction is the residual over the
      !! orbital-energy gaps, orthogonalised against the ones before it. One
      !! integral pass an iteration.
      type(czt_molecule_t), intent(in) :: mol
      type(spin_block_t), intent(in) :: sp(2)
      real(dp), intent(in) :: zero_h(:, :), bounds(:, :)
      real(dp), intent(in) :: x_a(:, :), x_b(:, :)
      real(dp), allocatable, intent(out) :: z_a(:, :), z_b(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gaps_a(:, :), gaps_b(:, :), rhs(:), resid(:), precond(:)
      real(dp), allocatable :: trial(:, :), image(:, :), small(:, :), proj(:), y(:)
      real(dp), allocatable :: u_a(:, :, :), u_b(:, :, :), au_a(:, :, :), au_b(:, :, :)
      real(dp) :: rnorm, norm
      integer :: n_a, n_b, n, k, m, pass
      logical :: singular
      character(len=LINE_LEN) :: line

      call gap_matrix(sp(1), gaps_a)
      call gap_matrix(sp(2), gaps_b)
      n_a = size(x_a)
      n_b = size(x_b)
      n = n_a + n_b
      allocate (z_a, mold=x_a)
      allocate (z_b, mold=x_b)
      z_a = 0.0_dp
      z_b = 0.0_dp
      if (n == 0) return

      allocate (rhs(n), resid(n), precond(n))
      rhs(1:n_a) = reshape(x_a, [n_a])
      rhs(n_a + 1:n) = reshape(x_b, [n_b])
      precond(1:n_a) = reshape(gaps_a, [n_a])
      precond(n_a + 1:n) = reshape(gaps_b, [n_b])
      allocate (trial(n, Z_MAX_ITER), image(n, Z_MAX_ITER))
      allocate (small(Z_MAX_ITER, Z_MAX_ITER), proj(Z_MAX_ITER))
      allocate (u_a(size(x_a, 1), size(x_a, 2), 1), u_b(size(x_b, 1), size(x_b, 2), 1))
      allocate (au_a, mold=u_a)
      allocate (au_b, mold=u_b)

      resid = rhs
      do k = 1, min(Z_MAX_ITER, n)
         rnorm = maxval(abs(resid))
         if (rnorm < Z_TOL) exit

         ! The next direction, orthogonalised twice: once loses digits when
         ! the residual is nearly inside the span already.
         trial(:, k) = resid/precond
         do pass = 1, 2
            do m = 1, k - 1
               trial(:, k) = trial(:, k) - dot_product(trial(:, m), trial(:, k))*trial(:, m)
            end do
         end do
         norm = sqrt(dot_product(trial(:, k), trial(:, k)))
         if (norm < tiny(1.0_dp)) then
            call error%set(ERROR_VALIDATION, "the unrestricted RI-MP2 gradient's "// &
                           "Z-vector solve ran out of directions before converging")
            return
         end if
         trial(:, k) = trial(:, k)/norm

         u_a(:, :, 1) = reshape(trial(1:n_a, k), shape(x_a))
         u_b(:, :, 1) = reshape(trial(n_a + 1:n, k), shape(x_b))
         call response_product_uhf(mol, sp(1)%c(:, 1:sp(1)%n_o), sp(1)%c(:, sp(1)%n_o + 1:), &
                                   sp(2)%c(:, 1:sp(2)%n_o), sp(2)%c(:, sp(2)%n_o + 1:), &
                                   gaps_a, gaps_b, zero_h, bounds, u_a, u_b, .false., &
                                   au_a, au_b, error)
         if (error%has_error()) return
         image(1:n_a, k) = reshape(au_a(:, :, 1), [n_a])
         image(n_a + 1:n, k) = reshape(au_b(:, :, 1), [n_b])

         do m = 1, k
            small(m, k) = dot_product(trial(:, m), image(:, k))
            small(k, m) = small(m, k)
         end do
         proj(k) = dot_product(trial(:, k), rhs)
         call solve_small(small(1:k, 1:k), proj(1:k), y, singular)
         if (singular) then
            call error%set(ERROR_VALIDATION, "the unrestricted RI-MP2 gradient's "// &
                           "Z-vector equations are singular: the reference has an "// &
                           "orbital rotation that costs no energy, and its response "// &
                           "is not defined")
            return
         end if
         resid = rhs - matmul(image(:, 1:k), y)
      end do
      rnorm = maxval(abs(resid))
      if (rnorm >= Z_TOL) then
         call error%set(ERROR_VALIDATION, "the unrestricted RI-MP2 gradient's Z-vector "// &
                        "equations did not converge")
         return
      end if
      if (allocated(y)) then
         resid = matmul(trial(:, 1:size(y)), y)
         z_a = reshape(resid(1:n_a), shape(x_a))
         z_b = reshape(resid(n_a + 1:n), shape(x_b))
         write (line, "(a,i0,a,es9.2)") "  unrestricted Z-vector: ", size(y), &
            " iterations, residual ", rnorm
         call logger%verbose(trim(line))
      end if
   end subroutine solve_z

   subroutine solve_small(matrix, rhs, x, singular)
      !! A small dense symmetric system by elimination with partial pivoting
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

   subroutine gap_matrix(spin, gaps)
      !! `eps_a - eps_i`, virtual by occupied
      type(spin_block_t), intent(in) :: spin
      real(dp), allocatable, intent(out) :: gaps(:, :)

      integer :: i, a

      allocate (gaps(spin%n_v, spin%n_o))
      do i = 1, spin%n_o
         do a = 1, spin%n_v
            gaps(a, i) = spin%eps(spin%n_o + a) - spin%eps(i)
         end do
      end do
   end subroutine gap_matrix

   subroutine weighted_density(spin, g, z, frozen)
      !! One spin's energy-weighted density, reference included
      !!
      !! The gradient contracts it with minus the overlap derivative, summed
      !! over every ordered pair of orbitals. Block by block, with `G` the
      !! reference's mean field of the relaxed density in this spin's
      !! orbitals:
      !!
      !!     correlated occupied pair  (l1_ij + l1_ji)/4 + (e_i + e_j) P_ij/2 + G_ij
      !!     any other occupied pair   G_ij, less z e_i/2 across frozen-correlated
      !!     virtual pair              (l1_ab + l1_ba)/4 + (e_a + e_b) P_ab/2
      !!     virtual-occupied          (l1_ia - z_ai e_i)/2
      !!
      !! and the reference's own `e_i` on each occupied diagonal.
      type(spin_block_t), intent(inout) :: spin
      real(dp), intent(in) :: g(:, :)   !! (n_ao, n_ao)
      real(dp), intent(in) :: z(:, :)   !! (n_v, n_o)
      integer, intent(in) :: frozen

      real(dp), allocatable :: g_mo(:, :)
      integer :: p, q, i, a, f, n_o, n_mo

      n_o = spin%n_o
      n_mo = size(spin%c, 2)
      spin%w = 0.0_dp
      if (n_o == 0) return
      g_mo = matmul(transpose(spin%c(:, 1:n_o)), matmul(g, spin%c(:, 1:n_o)))
      spin%w(1:n_o, 1:n_o) = g_mo
      do q = frozen + 1, n_o
         do p = frozen + 1, n_o
            spin%w(p, q) = spin%w(p, q) + 0.25_dp*(spin%l1(p, q) + spin%l1(q, p)) &
                           + 0.5_dp*(spin%eps(p) + spin%eps(q))*spin%p(p, q)
         end do
      end do
      ! `p(f, i)` holds minus half the frozen-correlated rotation, so this
      ! is `-z e_i/2`.
      do i = frozen + 1, n_o
         do f = 1, frozen
            spin%w(f, i) = spin%w(f, i) + spin%p(f, i)*spin%eps(i)
            spin%w(i, f) = spin%w(i, f) + spin%p(f, i)*spin%eps(i)
         end do
      end do
      do q = n_o + 1, n_mo
         do p = n_o + 1, n_mo
            spin%w(p, q) = 0.25_dp*(spin%l1(p, q) + spin%l1(q, p)) &
                           + 0.5_dp*(spin%eps(p) + spin%eps(q))*spin%p(p, q)
         end do
      end do
      do i = 1, n_o
         do a = 1, spin%n_v
            spin%w(n_o + a, i) = 0.5_dp*(spin%l1(i, n_o + a) - z(a, i)*spin%eps(i))
            spin%w(i, n_o + a) = spin%w(n_o + a, i)
         end do
      end do
      do i = 1, n_o
         spin%w(i, i) = spin%w(i, i) + spin%eps(i)
      end do
   end subroutine weighted_density

   subroutine fitted_densities(sp, jm12, n_ao, n_aux, gamma_ao, gamma_pq)
      !! The three- and two-index densities the fitted derivative contracts
      !!
      !! Summed over spin and halved, which is the closed-shell routine's
      !! convention: `fitted_two_particle_gradient` forms
      !! `4 Gamma (mu nu|P)^x - 2 gamma (P|Q)^x` and the functional here wants
      !! half of each.
      type(spin_block_t), intent(in) :: sp(2)
      real(dp), intent(in) :: jm12(:, :)
      integer, intent(in) :: n_ao, n_aux
      real(dp), allocatable, intent(out) :: gamma_ao(:, :, :), gamma_pq(:, :)

      real(dp), allocatable :: fit(:, :), gflat(:, :), c_act(:, :), c_vir(:, :)
      integer :: s, pp, i, a

      allocate (gamma_ao(n_ao, n_ao, n_aux), gamma_pq(n_aux, n_aux))
      gamma_ao = 0.0_dp
      gamma_pq = 0.0_dp
      do s = 1, 2
         if (sp(s)%n_ov == 0) cycle
         c_act = sp(s)%c(:, sp(s)%n_o - sp(s)%n_act + 1:sp(s)%n_o)
         c_vir = sp(s)%c(:, sp(s)%n_o + 1:)
         !$omp parallel do default(none) shared(gamma_ao, sp, s, c_act, c_vir, n_aux) &
         !$omp    private(pp)
         do pp = 1, n_aux
            gamma_ao(:, :, pp) = gamma_ao(:, :, pp) + 0.5_dp* &
                                 matmul(c_act, matmul(sp(s)%gamma(pp, :, :), transpose(c_vir)))
         end do
         !$omp end parallel do
         ! `gamma_PQ = sum_ia Gamma^P_ia c^Q_ia`, the fit coefficients being
         ! `B J^-1/2`.
         allocate (fit(sp(s)%n_ov, n_aux), gflat(sp(s)%n_ov, n_aux))
         call pic_gemm(sp(s)%b, jm12, fit, transb="T", beta=0.0_dp)
         do a = 1, sp(s)%n_v
            do i = 1, sp(s)%n_act
               gflat(i + (a - 1)*sp(s)%n_act, :) = sp(s)%gamma(:, i, a)
            end do
         end do
         call pic_gemm(gflat, fit, gamma_pq, transa="T", alpha=0.5_dp, beta=1.0_dp)
         deallocate (fit, gflat)
      end do
   end subroutine fitted_densities

end module mqc_czt_uri_mp2_gradient
