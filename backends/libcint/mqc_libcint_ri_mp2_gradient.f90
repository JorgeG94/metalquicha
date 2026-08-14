!! The RI-MP2 nuclear gradient over a closed-shell reference
module mqc_libcint_ri_mp2_gradient
   !! **What is differentiated.** An exact-ERI restricted Hartree-Fock
   !! reference with a *fitted* correlation energy -- which is what an `ri-mp2`
   !! deck asks for, the auxiliary basis naming the correlation fit and not the
   !! SCF. That is why the four-centre coupling tensor and the Fock derivative
   !! below are never fitted: approximating them would differentiate an energy
   !! nobody computed.
   !!
   !! The formulation is Stocks, Palethorpe and Barca, J. Chem. Theory Comput.
   !! 2024, 20, 2505, section 3. Only two things differ from the conventional
   !! MP2 gradient next door:
   !!
   !! * the two-particle density stays three-index, `Gamma^P_ia`, and contracts
   !!   against three- and two-centre derivative integrals rather than against
   !!   four-centre ones; and
   !! * the Lagrangian is built from `(pq|P)` instead of from a four-index
   !!   density.
   !!
   !! Everything one-particle -- the relaxed density, the energy-weighted
   !! density, the Z-vector, and the core-Hamiltonian, overlap and
   !! Hartree-Fock two-electron derivative terms -- is
   !! `mqc_libcint_mp2_gradient`'s, reused rather than rewritten.
   !!
   !! **Conventions that cost a day to establish**, none of them in the paper,
   !! all of them settled against the numpy prototype in
   !! `tools/cpu_validation/ri_mp2_gradient.py`:
   !!
   !! 1. The paper's `L'` and `L''` index the Lagrangian's own index first; the
   !!    assembly this feeds carries it second, so `imat` is their transpose.
   !!    The diagonals and the norms agree either way -- only the full matrix
   !!    comparison shows it.
   !! 2. The paper's `A_pqrs` contracted with a density is *twice* the response
   !!    operator the conventional assembly uses.
   !! 3. Terms 3 and 4 stayed wrong through three revisions while translational
   !!    invariance held at 1e-15. It is blind to every one of these.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, atom_ao_blocks, &
                                    build_df_mo_tensor, three_centre, two_centre, &
                                    metric_inverse_sqrt
   use mqc_libcint_cphf, only: cphf_solve
   use mqc_libcint_gradient, only: nuclear_repulsion_gradient, one_electron_deriv, &
                                   iprinv_deriv_at, three_centre_deriv, &
                                   two_centre_deriv, DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_libcint_direct, only: schwarz_bounds
   use mqc_libcint_mp2_gradient, only: gamma1_intermediates, two_electron_potential, &
                                       two_electron_mp2_terms, IN_CORE_LIMIT
   implicit none
   private

   public :: libcint_ri_mp2_gradient

contains

   subroutine libcint_ri_mp2_gradient(mol, aux, coeff, orbital_energies, n_occ, &
                                      gradient, error, n_frozen, force_direct)
      !! dE(RI-MP2)/dR for a closed-shell reference, in Hartree/Bohr
      type(libcint_molecule_t), intent(in) :: mol
      type(libcint_molecule_t), intent(in) :: aux   !! Auxiliary basis, same atoms
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: n_frozen
      logical, intent(in), optional :: force_direct
         !! Recompute the reference integrals rather than storing them, whatever
         !! the size. The check harness passes it so both paths are exercised on
         !! a case small enough to take either -- a threshold nothing ever
         !! crosses in a test is a threshold nothing ever tests.

      real(dp), allocatable :: bia(:, :, :), bia_raw(:, :, :), metric(:, :), jm12(:, :)
      real(dp), allocatable :: t2(:, :, :, :), xbar(:, :, :, :)
      real(dp), allocatable :: gamma(:, :, :), gamma_pq(:, :), gamma_ao(:, :, :)
      real(dp), allocatable :: doo(:, :), dvv(:, :), dm1mo(:, :), zeta(:, :)
      real(dp), allocatable :: imat(:, :), im1(:, :), im1_t(:, :)
      real(dp), allocatable :: lag_o(:, :), lag_v(:, :), pq_p(:, :, :)
      real(dp), allocatable :: c_occ(:, :), c_vir(:, :)
      real(dp), allocatable :: hf_density(:, :), dm1(:, :), dm1_total(:, :), dm1p(:, :)
      real(dp), allocatable :: veff(:, :), xvo(:, :, :), zvec(:, :, :)
      real(dp), allocatable :: work(:, :), work_t(:, :), p_occ(:, :), vhf_s1occ(:, :)
      real(dp), allocatable :: s1(:, :, :), h1(:, :, :), kin(:, :, :), zero_h(:, :)
      real(dp), allocatable :: vrinv(:, :, :), hcore_a(:, :, :)
      real(dp), allocatable :: de2(:, :), vhf1(:, :, :, :), eri(:, :, :, :)
      real(dp), allocatable :: bounds(:, :), gamma_null(:, :, :, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, n_mo, n_o, n_v, n_aux, frozen
      integer :: i, j, a, b, p, q, iatom, comp, p0, p1
      real(dp) :: denom
      logical :: dense

      n_ao = mol%nao
      n_mo = size(coeff, 2)
      n_o = n_occ
      n_v = n_mo - n_occ

      frozen = 0
      if (present(n_frozen)) frozen = n_frozen
      if (frozen /= 0) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 gradient does not implement a "// &
                        "frozen core: the relaxed density gains occupied-frozen and "// &
                        "virtual-frozen blocks, which are not built here.")
         return
      end if
      if (n_v < 1 .or. n_o < 1) then
         call error%set(ERROR_VALIDATION, "the RI-MP2 gradient needs both occupied "// &
                        "and virtual orbitals")
         return
      end if

      ! Stored or recomputed, for the *reference* integrals only -- the ones the
      ! Z-vector's operator and the reference potential are built from. Those
      ! stay exact whatever the correlation is fitted with, because an `ri-mp2`
      ! energy fits nothing in its SCF, so the response has to be the exact
      ! reference's or it answers a question nobody asked.
      !
      ! What is fitted is not what makes this decision. Storing `n_ao^4` doubles
      ! would put the ceiling of a method chosen to avoid `n_ao^4` back where RI
      ! removed it, so past the limit the same integrals are recomputed instead
      ! and nothing about the answer changes.
      dense = real(n_ao, dp)**4*8.0_dp <= IN_CORE_LIMIT
      if (present(force_direct)) then
         if (force_direct) dense = .false.
      end if

      allocate (c_occ(n_ao, n_o), c_vir(n_ao, n_v))
      c_occ = coeff(:, 1:n_o)
      c_vir = coeff(:, n_o + 1:n_mo)

      ! ---- the fitted integrals, and the amplitudes over them --------------
      ! `build_df_mo_tensor` lays its result out `(a, P, i)`, the shape the
      ! energy's gemms want. Everything below indexes `(i, a, P)` instead,
      ! matching the equations and the numpy reference, so it is repacked once
      ! here rather than transposed at every use.
      call build_df_mo_tensor(mol, aux, c_occ, c_vir, bia_raw, error)
      if (error%has_error()) return
      n_aux = size(bia_raw, 2)
      allocate (bia(n_o, n_v, n_aux))
      do p = 1, n_aux
         do i = 1, n_o
            bia(i, :, p) = bia_raw(:, p, i)
         end do
      end do
      deallocate (bia_raw)

      call two_centre(aux, metric)
      call metric_inverse_sqrt(metric, jm12, error)
      if (error%has_error()) return

      allocate (t2(n_o, n_o, n_v, n_v))
      do b = 1, n_v
         do a = 1, n_v
            do j = 1, n_o
               do i = 1, n_o
                  denom = orbital_energies(i) + orbital_energies(j) &
                          - orbital_energies(n_occ + a) - orbital_energies(n_occ + b)
                  t2(i, j, a, b) = sum(bia(i, a, :)*bia(j, b, :))/denom
               end do
            end do
         end do
      end do

      ! ---- the three- and two-index densities (eqs 8, 10) ------------------
      allocate (xbar(n_o, n_o, n_v, n_v))
      xbar = 2.0_dp*t2
      do b = 1, n_v
         do a = 1, n_v
            xbar(:, :, a, b) = xbar(:, :, a, b) - t2(:, :, b, a)
         end do
      end do

      call build_gamma(xbar, bia, jm12, n_o, n_v, n_aux, gamma)
      call build_gamma_metric(gamma, bia, jm12, n_o, n_v, n_aux, gamma_pq)
      deallocate (xbar)

      ! ---- the unrelaxed one-particle blocks -------------------------------
      call gamma1_intermediates(t2, n_o, n_v, doo, dvv)
      deallocate (t2)

      allocate (dm1mo(n_mo, n_mo))
      dm1mo = 0.0_dp
      dm1mo(1:n_o, 1:n_o) = doo + transpose(doo)
      dm1mo(n_o + 1:n_mo, n_o + 1:n_mo) = dvv + transpose(dvv)

      ! ---- the Lagrangian (eqs 13, 14) -------------------------------------
      !
      ! `(pq|P)` over every molecular orbital, so the three-centre integrals are
      ! transformed in full rather than only to the occupied-virtual block.
      call three_centre(mol, aux, work)
      allocate (pq_p(n_mo, n_mo, n_aux))
      call transform_three_centre(work, coeff, n_ao, n_mo, n_aux, pq_p)
      deallocate (work)

      allocate (lag_o(n_o, n_mo), lag_v(n_v, n_mo))
      lag_o = 0.0_dp
      lag_v = 0.0_dp
      do q = 1, n_mo
         do a = 1, n_v
            do i = 1, n_o
               lag_o(i, q) = lag_o(i, q) &
                             + 2.0_dp*sum(gamma(:, i, a)*pq_p(q, n_o + a, :))
            end do
         end do
         do a = 1, n_v
            do i = 1, n_o
               lag_v(a, q) = lag_v(a, q) &
                             + 2.0_dp*sum(gamma(:, i, a)*pq_p(i, q, :))
            end do
         end do
      end do
      deallocate (pq_p)

      ! Transposed into the slot: the paper indexes the Lagrangian's own index
      ! first, the assembly below carries it second.
      allocate (imat(n_mo, n_mo))
      imat(:, 1:n_o) = -transpose(lag_o)
      imat(:, n_o + 1:n_mo) = -transpose(lag_v)

      ! ---- the Z-vector ----------------------------------------------------
      allocate (hf_density(n_ao, n_ao))
      hf_density = 2.0_dp*matmul(c_occ, transpose(c_occ))

      if (dense) then
         call mol%eris(eri)
         allocate (bounds(0, 0))
      else
         allocate (eri(0, 0, 0, 0))
         call schwarz_bounds(mol, bounds, error)
         if (error%has_error()) return
      end if
      allocate (zero_h(n_ao, n_ao), veff(n_ao, n_ao))
      zero_h = 0.0_dp

      allocate (dm1(n_ao, n_ao))
      dm1 = matmul(coeff, matmul(dm1mo, transpose(coeff)))
      call two_electron_potential(dense, mol, eri, bounds, zero_h, dm1, veff, error)
      if (error%has_error()) return
      veff = 2.0_dp*veff

      allocate (xvo(n_v, n_o, 1))
      xvo(:, :, 1) = matmul(transpose(c_vir), matmul(veff, c_occ))
      do i = 1, n_o
         do a = 1, n_v
            xvo(a, i, 1) = xvo(a, i, 1) + imat(i, n_o + a) - imat(n_o + a, i)
         end do
      end do

      ! Tighter than the solver's default, for the reason the conventional
      ! gradient gives: the Z-vector residual enters the answer directly rather
      ! than quadratically. Handing over the tensor when there is one saves a
      ! second identical pass; when there is not, the solver goes direct with
      ! the same Schwarz bound the potentials above used.
      if (dense) then
         call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zvec, &
                         error=error, mo_rhs=xvo, tol=1.0e-12_dp, max_iter=200, &
                         eri_in=eri)
      else
         call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zvec, &
                         error=error, mo_rhs=xvo, tol=1.0e-12_dp, max_iter=200)
      end if
      if (error%has_error()) return

      dm1mo(n_o + 1:n_mo, 1:n_o) = dm1mo(n_o + 1:n_mo, 1:n_o) + zvec(:, :, 1)
      dm1mo(1:n_o, n_o + 1:n_mo) = dm1mo(1:n_o, n_o + 1:n_mo) + transpose(zvec(:, :, 1))

      imat(n_o + 1:n_mo, 1:n_o) = transpose(imat(1:n_o, n_o + 1:n_mo))
      allocate (im1(n_ao, n_ao))
      im1 = matmul(coeff, matmul(imat, transpose(coeff)))

      ! ---- the energy-weighted density -------------------------------------
      allocate (zeta(n_mo, n_mo))
      do q = 1, n_mo
         do p = 1, n_mo
            zeta(p, q) = 0.5_dp*(orbital_energies(p) + orbital_energies(q))
         end do
      end do
      do i = 1, n_o
         do a = 1, n_v
            zeta(n_o + a, i) = orbital_energies(i)
            zeta(i, n_o + a) = orbital_energies(i)
         end do
      end do
      zeta = zeta*dm1mo
      allocate (work(n_ao, n_ao))
      work = matmul(coeff, matmul(zeta, transpose(coeff)))
      deallocate (zeta)

      dm1 = matmul(coeff, matmul(dm1mo, transpose(coeff)))
      allocate (p_occ(n_ao, n_ao))
      p_occ = matmul(c_occ, transpose(c_occ))
      call two_electron_potential(dense, mol, eri, bounds, zero_h, &
                                  dm1 + transpose(dm1), veff, error)
      if (error%has_error()) return
      allocate (vhf_s1occ(n_ao, n_ao))
      vhf_s1occ = matmul(p_occ, matmul(veff, p_occ))

      allocate (dm1p(n_ao, n_ao), dm1_total(n_ao, n_ao))
      dm1p = hf_density + 2.0_dp*dm1
      dm1_total = hf_density + dm1

      do i = 1, n_o
         do q = 1, n_ao
            do p = 1, n_ao
               work(p, q) = work(p, q) + 2.0_dp*orbital_energies(i)*c_occ(p, i)*c_occ(q, i)
            end do
         end do
      end do

      ! ---- assembly --------------------------------------------------------
      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp
      call nuclear_repulsion_gradient(mol, gradient)

      allocate (de2(3, mol%natm), vhf1(n_ao, n_ao, 3, mol%natm))
      de2 = 0.0_dp
      vhf1 = 0.0_dp
      ! Only the reference-density matrices: the fitted two-particle density
      ! contracts against three-centre derivatives instead, below.
      allocate (gamma_null(0, 0, 0, 0))
      call two_electron_mp2_terms(mol, gamma_null, hf_density, de2, vhf1, &
                                  with_gamma=.false.)

      call build_gamma_ao(gamma, c_occ, c_vir, n_ao, n_o, n_v, n_aux, gamma_ao)
      call fitted_two_particle_gradient(mol, aux, gamma_ao, gamma_pq, gradient)

      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (vrinv(n_ao, n_ao, 3), hcore_a(n_ao, n_ao, 3))
      allocate (im1_t(n_ao, n_ao), work_t(n_ao, n_ao))
      im1_t = transpose(im1)
      work_t = transpose(work)

      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         call iprinv_deriv_at(mol, iatom, vrinv)
         vrinv = -mol%charges(iatom)*vrinv
         if (counts(iatom) > 0) then
            vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         end if
         do comp = 1, 3
            hcore_a(:, :, comp) = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(hcore_a(:, :, comp)*dm1_total) &
                                    - sum(vhf1(:, :, comp, iatom)*dm1p)
         end do

         if (counts(iatom) == 0) cycle
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(s1(p0:p1, :, comp)*im1(p0:p1, :)) &
                                    + sum(s1(p0:p1, :, comp)*im1_t(p0:p1, :)) &
                                    - sum(s1(p0:p1, :, comp)*work(p0:p1, :)) &
                                    - sum(s1(p0:p1, :, comp)*work_t(p0:p1, :)) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*vhf_s1occ(p0:p1, :))
         end do
      end do

   end subroutine libcint_ri_mp2_gradient

   subroutine build_gamma(xbar, bia, jm12, n_o, n_v, n_aux, gamma)
      !! `Gamma^P_ia = sum_{bjQ} (2 X^ab_ij - X^ba_ij) B^Q_jb J^{-1/2}_PQ` (eq 8)
      real(dp), intent(in) :: xbar(:, :, :, :), bia(:, :, :), jm12(:, :)
      integer, intent(in) :: n_o, n_v, n_aux
      real(dp), allocatable, intent(out) :: gamma(:, :, :)

      real(dp), allocatable :: half(:, :, :)
      integer :: i, j, a, b, qq

      allocate (half(n_aux, n_o, n_v), gamma(n_aux, n_o, n_v))
      half = 0.0_dp
      !$omp parallel do collapse(2) default(none) &
      !$omp    shared(half, xbar, bia, n_o, n_v, n_aux) private(i, a, j, b, qq)
      do a = 1, n_v
         do i = 1, n_o
            do b = 1, n_v
               do j = 1, n_o
                  do qq = 1, n_aux
                     half(qq, i, a) = half(qq, i, a) + xbar(i, j, a, b)*bia(j, b, qq)
                  end do
               end do
            end do
         end do
      end do
      !$omp end parallel do

      do a = 1, n_v
         do i = 1, n_o
            gamma(:, i, a) = matmul(jm12, half(:, i, a))
         end do
      end do
      deallocate (half)
   end subroutine build_gamma

   subroutine build_gamma_metric(gamma, bia, jm12, n_o, n_v, n_aux, gamma_pq)
      !! `gamma_PQ = sum_{iaR} Gamma^P_ia B^R_ia J^{-1/2}_QR` (eq 10)
      real(dp), intent(in) :: gamma(:, :, :), bia(:, :, :), jm12(:, :)
      integer, intent(in) :: n_o, n_v, n_aux
      real(dp), allocatable, intent(out) :: gamma_pq(:, :)

      real(dp), allocatable :: half(:, :)
      integer :: i, a

      allocate (half(n_aux, n_aux), gamma_pq(n_aux, n_aux))
      half = 0.0_dp
      do a = 1, n_v
         do i = 1, n_o
            call outer_accumulate(gamma(:, i, a), bia(i, a, :), half)
         end do
      end do
      gamma_pq = matmul(half, transpose(jm12))
      deallocate (half)
   end subroutine build_gamma_metric

   subroutine outer_accumulate(u, v, m)
      !! `m = m + u v^T`
      real(dp), intent(in) :: u(:), v(:)
      real(dp), intent(inout) :: m(:, :)
      integer :: k

      do k = 1, size(v)
         m(:, k) = m(:, k) + u*v(k)
      end do
   end subroutine outer_accumulate

   subroutine build_gamma_ao(gamma, c_occ, c_vir, n_ao, n_o, n_v, n_aux, gamma_ao)
      !! `Gamma^P_munu`, the three-index density back-transformed
      real(dp), intent(in) :: gamma(:, :, :), c_occ(:, :), c_vir(:, :)
      integer, intent(in) :: n_ao, n_o, n_v, n_aux
      real(dp), allocatable, intent(out) :: gamma_ao(:, :, :)

      integer :: pp

      allocate (gamma_ao(n_ao, n_ao, n_aux))
      !$omp parallel do default(none) &
      !$omp    shared(gamma_ao, gamma, c_occ, c_vir, n_aux) private(pp)
      do pp = 1, n_aux
         gamma_ao(:, :, pp) = matmul(c_occ, matmul(gamma(pp, :, :), transpose(c_vir)))
      end do
      !$omp end parallel do
   end subroutine build_gamma_ao

   subroutine transform_three_centre(three, coeff, n_ao, n_mo, n_aux, pq_p)
      !! `(mu nu|P)` to `(pq|P)`, both molecular indices
      real(dp), intent(in) :: three(:, :), coeff(:, :)
      integer, intent(in) :: n_ao, n_mo, n_aux
      real(dp), intent(out) :: pq_p(:, :, :)

      real(dp), allocatable :: block(:, :)
      integer :: pp

      allocate (block(n_ao, n_ao))
      !$omp parallel do default(none) &
      !$omp    shared(pq_p, three, coeff, n_ao, n_aux) private(pp, block)
      do pp = 1, n_aux
         block = reshape(three(:, pp), [n_ao, n_ao])
         pq_p(:, :, pp) = matmul(transpose(coeff), matmul(block, coeff))
      end do
      !$omp end parallel do
      deallocate (block)
   end subroutine transform_three_centre

   subroutine fitted_two_particle_gradient(mol, aux, gamma_ao, gamma_pq, gradient)
      !! The first two terms of eq 6, the only ones the fitting changes
      !!
      !! `4 sum Gamma^P_munu (P|munu)^xi` over both the orbital centres and the
      !! auxiliary one, then `-2 sum gamma_PQ (P|Q)^xi`. libcint's `ip`
      !! integrals differentiate the electronic coordinate, so the derivative
      !! with respect to a nucleus carrying that function is minus them.
      type(libcint_molecule_t), intent(in) :: mol, aux
      real(dp), intent(in) :: gamma_ao(:, :, :), gamma_pq(:, :)
      real(dp), intent(inout) :: gradient(:, :)

      real(dp), allocatable :: ip1(:, :, :, :), ip2(:, :, :, :), j1(:, :, :)
      integer, allocatable :: offsets(:), counts(:), aux_offsets(:), aux_counts(:)
      integer :: iatom, comp, p0, p1, q0, q1

      allocate (offsets(mol%natm), counts(mol%natm))
      allocate (aux_offsets(aux%natm), aux_counts(aux%natm))
      call atom_ao_blocks(mol, offsets, counts)
      call atom_ao_blocks(aux, aux_offsets, aux_counts)

      call three_centre_deriv(mol, aux, 1, ip1)
      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)
         if (counts(iatom) == 0) cycle
         do comp = 1, 3
            ! The bra index on this atom, and then the ket: `(mu nu|P)` is
            ! symmetric in mu and nu, so the second is the first transposed.
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - 4.0_dp*sum(ip1(p0:p1, :, :, comp) &
                                                 *gamma_ao(p0:p1, :, :)) &
                                    - 4.0_dp*sum(ip1(p0:p1, :, :, comp) &
                                                 *transpose_first_two(gamma_ao, p0, p1))
         end do
      end do
      deallocate (ip1)

      call three_centre_deriv(mol, aux, 2, ip2)
      do iatom = 1, aux%natm
         q0 = aux_offsets(iatom) + 1
         q1 = aux_offsets(iatom) + aux_counts(iatom)
         if (aux_counts(iatom) == 0) cycle
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    - 4.0_dp*sum(ip2(:, :, q0:q1, comp) &
                                                 *gamma_ao(:, :, q0:q1))
         end do
      end do
      deallocate (ip2)

      call two_centre_deriv(aux, j1)
      do iatom = 1, aux%natm
         q0 = aux_offsets(iatom) + 1
         q1 = aux_offsets(iatom) + aux_counts(iatom)
         if (aux_counts(iatom) == 0) cycle
         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + 2.0_dp*sum(j1(q0:q1, :, comp)*gamma_pq(q0:q1, :)) &
                                    + 2.0_dp*sum(j1(q0:q1, :, comp) &
                                                 *transpose(gamma_pq(:, q0:q1)))
         end do
      end do
      deallocate (j1)
   end subroutine fitted_two_particle_gradient

   function transpose_first_two(g, p0, p1) result(t)
      !! `g(nu, mu, P)` for `mu` in a block, as the ket-side contraction needs
      real(dp), intent(in) :: g(:, :, :)
      integer, intent(in) :: p0, p1
      real(dp), allocatable :: t(:, :, :)
      integer :: pp

      allocate (t(p1 - p0 + 1, size(g, 1), size(g, 3)))
      do pp = 1, size(g, 3)
         t(:, :, pp) = transpose(g(:, p0:p1, pp))
      end do
   end function transpose_first_two

end module mqc_libcint_ri_mp2_gradient
