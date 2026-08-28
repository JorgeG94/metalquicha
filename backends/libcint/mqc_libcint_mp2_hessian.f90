!! The fixed-density second-derivative skeleton of the MP2 correlation Hessian
module mqc_libcint_mp2_hessian
   !! Unit 1.3 of the MP2 Hessian ladder (`mp2-hessian-phased-plan.md`): the
   !! per-atom-pair scalar
   !!
   !!     s = Gamma_eff . (mu nu|la si)^{XY} + D_rel . h^{XY} + I . S^{XY}
   !!
   !! contracted against the **unperturbed** densities -- pass 1 of the
   !! two-pass assembly, mirrored from pycc's `_hessian_blocks` (`route ==
   !! 'aod'`). The orbital response and the perturbed densities are later
   !! units and deliberately absent here.
   !!
   !! **The reference skeleton rides in the same sweep.** The second-derivative
   !! two-electron integrals dominate the cost of a correlated Hessian, and the
   !! SCF reference needs the same ones contracted against its own separable
   !! density. Generating them once and depositing into two accumulators is
   !! pycc's `Href` fold, and it is also what makes the reference block
   !! checkable against `partial_hessian`, which walks its own quartets.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, atom_ao_blocks
   use mqc_libcint_hess_ints, only: hess_1e_block, hess_rinv_block, hess_2e_block, &
                                    HESS_KIN_II, HESS_KIN_IJ, HESS_NUC_II, HESS_NUC_IJ, &
                                    HESS_OVLP_II, HESS_OVLP_IJ, HESS_RINV_II, HESS_RINV_IJ, &
                                    HESS_ERI_II, HESS_ERI_IJ, HESS_ERI_IK
   implicit none
   private

   public :: mp2_skeleton_hessian

contains

   subroutine mp2_skeleton_hessian(mol, gamma_eff_ao, relaxed_density_mo, &
                                   energy_weighted_ao, coeff, orbital_energies, &
                                   n_occ, hess_corr, hess_ref, error)
      !! Both fixed-density skeletons, `(3, 3, natm, natm)` each
      !!
      !! `hess_corr(a, b, A, B)` is pycc's pass-1 scalar `s` for the atom pair
      !! `(A, B)` and Cartesian components `(a, b)`; `hess_ref` is the SCF
      !! skeleton over the same integrals -- the `partial_hessian` numbers,
      !! deposited from this sweep rather than a second one.
      !!
      !! **`energy_weighted_ao` is the total, and the trap is documented.**
      !! The gradient's `W` carries the reference's energy-weighted density,
      !! which is most of its norm; pycc's `I` is correlation only. The split
      !! is done here, once, rather than left to every caller:
      !!
      !!     I_corr = (W_total + 2 W_ref) / 2,   W_ref = 2 sum_i eps_i C_i C_i^T
      !!
      !! measured in `tools/mp2_hessian_oracle/CONVENTIONS.md` s.4a. Handing
      !! this routine a correlation-only `W` would double-count the reference's
      !! share -- give it exactly what `libcint_mp2_gradient` returns.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: gamma_eff_ao(:, :, :, :)
         !! From `build_effective_2pdm_ao`: chemist ordered, bra-ket
         !! symmetric, ket pair in libcint's natural order, with the
         !! two-electron half of `D_rel f^{XY}` already folded in.
      real(dp), intent(in) :: relaxed_density_mo(:, :)
         !! The relaxed one-particle density from the gradient, MO basis.
      real(dp), intent(in) :: energy_weighted_ao(:, :)
         !! The **total** energy-weighted density from the gradient -- the
         !! matrix that multiplies `dS/dR` there, reference share included.
      real(dp), intent(in) :: coeff(:, :)            !! C, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)    !! (n_mo), Hartree
      integer, intent(in) :: n_occ                   !! Doubly occupied count
      real(dp), allocatable, intent(out) :: hess_corr(:, :, :, :)
      real(dp), allocatable, intent(out) :: hess_ref(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: drel_ao(:, :), wcorr(:, :), wref(:, :), dref(:, :)
      real(dp), allocatable :: h2(:, :, :), hab(:, :, :), s2(:, :, :), sab(:, :, :)
      real(dp), allocatable :: tmp(:, :, :), r2(:, :, :), rab(:, :, :)
      real(dp), allocatable :: cross_c(:, :, :), cross_r(:, :, :)
      real(dp), allocatable :: eri_ii(:, :, :, :, :), eri_ij(:, :, :, :, :)
      real(dp), allocatable :: eri_ik(:, :, :, :, :)
      integer, allocatable :: owner(:), offsets(:), counts(:)
      real(dp) :: total_c(3, 3), total_r(3, 3)
      real(dp) :: wc, wr, w4c, w8c, gref, zc
      integer :: n_ao, n_mo, natm, iao, jao, kao, lao, ia, ja, ka, la
      integer :: a, b, comp, c, i

      if (error%has_error()) return

      n_ao = mol%nao
      n_mo = size(coeff, 2)
      natm = mol%natm
      if (size(gamma_eff_ao, 1) /= n_ao .or. size(relaxed_density_mo, 1) /= n_mo &
          .or. size(energy_weighted_ao, 1) /= n_ao &
          .or. size(orbital_energies) /= n_mo) then
         call error%set(ERROR_VALIDATION, "the MP2 skeleton Hessian's densities "// &
                        "and the basis disagree about their sizes.")
         return
      end if

      ! ---- the four densities the sweep contracts ---------------------------
      allocate (drel_ao(n_ao, n_ao))
      drel_ao = matmul(coeff, matmul(relaxed_density_mo, transpose(coeff)))

      allocate (wref(n_ao, n_ao))
      wref = 0.0_dp
      do i = 1, n_occ
         wref = wref + 2.0_dp*orbital_energies(i) &
                *matmul(reshape(coeff(:, i), [n_ao, 1]), &
                        reshape(coeff(:, i), [1, n_ao]))
      end do
      allocate (wcorr(n_ao, n_ao))
      wcorr = 0.5_dp*(energy_weighted_ao + 2.0_dp*wref)

      allocate (dref(n_ao, n_ao))
      dref = 2.0_dp*matmul(coeff(:, 1:n_occ), transpose(coeff(:, 1:n_occ)))

      allocate (offsets(natm), counts(natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (owner(n_ao))
      owner = 0
      do c = 1, natm
         do iao = offsets(c) + 1, offsets(c) + counts(c)
            owner(iao) = c
         end do
      end do

      allocate (hess_corr(3, 3, natm, natm), hess_ref(3, 3, natm, natm))
      hess_corr = 0.0_dp
      hess_ref = 0.0_dp

      ! ---- one electron, both derivatives on basis centres ------------------
      !
      ! `partial_hessian`'s deposits with two sets of densities instead of one.
      ! The factor of two is bra and ket contributing the same number; the
      ! signs follow the energy expression -- `+ D h^{XY}` against
      ! `- W_ref S^{XY}` for the reference, while the correlation's
      ! energy-weighted term enters with pycc's sign, `+ I S^{XY}`.
      call hess_1e_block(mol, HESS_KIN_II, h2, error)
      call hess_1e_block(mol, HESS_NUC_II, tmp, error)
      if (error%has_error()) return
      h2 = h2 + tmp
      deallocate (tmp)

      call hess_1e_block(mol, HESS_KIN_IJ, hab, error)
      call hess_1e_block(mol, HESS_NUC_IJ, tmp, error)
      if (error%has_error()) return
      hab = hab + tmp
      deallocate (tmp)

      call hess_1e_block(mol, HESS_OVLP_II, s2, error)
      call hess_1e_block(mol, HESS_OVLP_IJ, sab, error)
      if (error%has_error()) return

      ! **Statement for statement `partial_hessian`'s deposits**, not merely
      ! term for term. The second-derivative one-electron integrals on a steep
      ! core exponent reach ~1e5, so the accumulator swings through values that
      ! large before cancelling to order one -- and any regrouping (one fused
      ! update instead of an add and a subtract) moves the reference block by
      ! ~1e-11 absolute, which is exactly the divergence the 1e-13 gate against
      ! `partial_hessian` exists to catch. Measured here before the statements
      ! were split to match; do not fuse them back for tidiness.
      do comp = 1, 9
         a = (comp - 1)/3 + 1
         b = comp - 3*(a - 1)
         do jao = 1, n_ao
            ja = owner(jao)
            do iao = 1, n_ao
               ia = owner(iao)
               wc = drel_ao(iao, jao)
               hess_corr(a, b, ia, ia) = hess_corr(a, b, ia, ia) + 2.0_dp*wc*h2(iao, jao, comp)
               hess_corr(a, b, ia, ja) = hess_corr(a, b, ia, ja) + 2.0_dp*wc*hab(iao, jao, comp)
               wc = wcorr(iao, jao)
               hess_corr(a, b, ia, ia) = hess_corr(a, b, ia, ia) + 2.0_dp*wc*s2(iao, jao, comp)
               hess_corr(a, b, ia, ja) = hess_corr(a, b, ia, ja) + 2.0_dp*wc*sab(iao, jao, comp)
               wr = dref(iao, jao)
               hess_ref(a, b, ia, ia) = hess_ref(a, b, ia, ia) + 2.0_dp*wr*h2(iao, jao, comp)
               hess_ref(a, b, ia, ja) = hess_ref(a, b, ia, ja) + 2.0_dp*wr*hab(iao, jao, comp)
               wr = wref(iao, jao)
               hess_ref(a, b, ia, ia) = hess_ref(a, b, ia, ia) - 2.0_dp*wr*s2(iao, jao, comp)
               hess_ref(a, b, ia, ja) = hess_ref(a, b, ia, ja) - 2.0_dp*wr*sab(iao, jao, comp)
            end do
         end do
      end do
      deallocate (h2, hab, s2, sab)

      ! ---- one electron, at least one derivative on a nucleus ---------------
      !
      ! The operator-centre terms of the nuclear attraction, exactly as
      ! `partial_hessian` assembles them, run once per nucleus with both
      ! densities deposited from the same integral blocks.
      allocate (cross_c(3, 3, natm), cross_r(3, 3, natm))
      do c = 1, natm
         call hess_rinv_block(mol, c, HESS_RINV_II, r2, error)
         call hess_rinv_block(mol, c, HESS_RINV_IJ, rab, error)
         if (error%has_error()) return
         zc = mol%charges(c)

         cross_c = 0.0_dp
         cross_r = 0.0_dp
         do comp = 1, 9
            a = (comp - 1)/3 + 1
            b = comp - 3*(a - 1)
            do jao = 1, n_ao
               do iao = 1, n_ao
                  ia = owner(iao)
                  wc = r2(iao, jao, comp) + rab(iao, jao, comp)
                  cross_c(a, b, ia) = cross_c(a, b, ia) + drel_ao(iao, jao)*wc
                  cross_r(a, b, ia) = cross_r(a, b, ia) + dref(iao, jao)*wc
               end do
            end do
         end do
         deallocate (r2, rab)

         ! `-Z_C` is the charge and sign of the electron-nucleus attraction,
         ! which `hess_rinv_block` deliberately leaves off.
         total_c = 0.0_dp
         total_r = 0.0_dp
         do ia = 1, natm
            total_c = total_c + cross_c(:, :, ia)
            total_r = total_r + cross_r(:, :, ia)
            do b = 1, 3
               do a = 1, 3
                  hess_corr(a, b, ia, c) = hess_corr(a, b, ia, c) &
                                           + 2.0_dp*zc*cross_c(a, b, ia)
                  hess_corr(a, b, c, ia) = hess_corr(a, b, c, ia) &
                                           + 2.0_dp*zc*cross_c(b, a, ia)
                  hess_ref(a, b, ia, c) = hess_ref(a, b, ia, c) &
                                          + 2.0_dp*zc*cross_r(a, b, ia)
                  hess_ref(a, b, c, ia) = hess_ref(a, b, c, ia) &
                                          + 2.0_dp*zc*cross_r(b, a, ia)
               end do
            end do
         end do
         ! Both derivatives on the same nucleus: translational invariance of
         ! `1/|r - R_C|` in its three arguments folds it onto the sum of the
         ! basis-derivative terms.
         hess_corr(:, :, c, c) = hess_corr(:, :, c, c) - 2.0_dp*zc*total_c
         hess_ref(:, :, c, c) = hess_ref(:, :, c, c) - 2.0_dp*zc*total_r
      end do
      deallocate (cross_c, cross_r)

      ! ---- two electron -----------------------------------------------------
      !
      ! The sixteen ordered derivative placements of `d2(mn|ls)/dX dY` fold
      ! onto the three integrals libcint provides, once the loop runs over
      ! every quartet in every order. For a weight with the integral's full
      ! eightfold symmetry the folded coefficients are the familiar 4, 4 and 8;
      ! `Gamma_eff` only has bra-ket symmetry, so the general fold is written
      ! out instead: the four orbit terms below multiply `ipip1` and `ipvip1`,
      ! all eight multiply `ip1ip2`. Setting every orbit term equal recovers
      ! 4/4/8, which is the check that the fold is the same one
      ! `hess_2e_contract` uses.
      !
      ! Materialised whole for now: at the sizes the ladder gates on, the three
      ! arrays are a few megabytes, and the quartet-driven version that
      ! evaluates and forgets is the follow-up commit, not this one.
      call hess_2e_block(mol, HESS_ERI_II, eri_ii, error)
      call hess_2e_block(mol, HESS_ERI_IJ, eri_ij, error)
      call hess_2e_block(mol, HESS_ERI_IK, eri_ik, error)
      if (error%has_error()) return

      do lao = 1, n_ao
         la = owner(lao)
         do kao = 1, n_ao
            ka = owner(kao)
            do jao = 1, n_ao
               ja = owner(jao)
               do iao = 1, n_ao
                  ia = owner(iao)
                  w4c = gamma_eff_ao(iao, jao, kao, lao) &
                        + gamma_eff_ao(jao, iao, kao, lao) &
                        + gamma_eff_ao(kao, lao, iao, jao) &
                        + gamma_eff_ao(kao, lao, jao, iao)
                  w8c = w4c + gamma_eff_ao(iao, jao, lao, kao) &
                        + gamma_eff_ao(jao, iao, lao, kao) &
                        + gamma_eff_ao(lao, kao, iao, jao) &
                        + gamma_eff_ao(lao, kao, jao, iao)
                  gref = 0.5_dp*dref(iao, jao)*dref(kao, lao) &
                         - 0.125_dp*(dref(iao, lao)*dref(kao, jao) &
                                     + dref(iao, kao)*dref(jao, lao))
                  do comp = 1, 9
                     a = (comp - 1)/3 + 1
                     b = comp - 3*(a - 1)
                     hess_corr(a, b, ia, ia) = hess_corr(a, b, ia, ia) &
                                               + w4c*eri_ii(iao, jao, kao, lao, comp)
                     hess_corr(a, b, ia, ja) = hess_corr(a, b, ia, ja) &
                                               + w4c*eri_ij(iao, jao, kao, lao, comp)
                     hess_corr(a, b, ia, ka) = hess_corr(a, b, ia, ka) &
                                               + w8c*eri_ik(iao, jao, kao, lao, comp)
                     hess_ref(a, b, ia, ia) = hess_ref(a, b, ia, ia) &
                                              + 4.0_dp*gref*eri_ii(iao, jao, kao, lao, comp)
                     hess_ref(a, b, ia, ja) = hess_ref(a, b, ia, ja) &
                                              + 4.0_dp*gref*eri_ij(iao, jao, kao, lao, comp)
                     hess_ref(a, b, ia, ka) = hess_ref(a, b, ia, ka) &
                                              + 8.0_dp*gref*eri_ik(iao, jao, kao, lao, comp)
                  end do
               end do
            end do
         end do
      end do
      deallocate (eri_ii, eri_ij, eri_ik)
      deallocate (drel_ao, wcorr, wref, dref, owner, offsets, counts)
   end subroutine mp2_skeleton_hessian

end module mqc_libcint_mp2_hessian
