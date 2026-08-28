!! The fixed-density second-derivative skeleton of the MP2 correlation Hessian
module mqc_libcint_mp2_hessian
   !! Units 1.3-1.6 of the MP2 Hessian ladder (`mp2-hessian-phased-plan.md`).
   !!
   !! Unit 1.3, `mp2_skeleton_hessian`: the per-atom-pair scalar
   !!
   !!     s = Gamma_eff . (mu nu|la si)^{XY} + D_rel . h^{XY} + I . S^{XY}
   !!
   !! contracted against the **unperturbed** densities -- pass 1 of the
   !! two-pass assembly, mirrored from pycc's `_hessian_blocks` (`route ==
   !! 'aod'`). The perturbed densities are later units and deliberately
   !! absent here.
   !!
   !! **The reference skeleton rides in the same sweep.** The second-derivative
   !! two-electron integrals dominate the cost of a correlated Hessian, and the
   !! SCF reference needs the same ones contracted against its own separable
   !! density. Generating them once and depositing into two accumulators is
   !! pycc's `Href` fold, and it is also what makes the reference block
   !! checkable against `partial_hessian`, which walks its own quartets.
   !!
   !! Units 1.4-1.6: the per-perturbation first-order skeletons `f^(X)`,
   !! `S^(X)`, `<pq|rs>^(X)`; the skeleton Lagrangian `I'^(X)` and its
   !! orbital-response carriers; and the orbital-response term of pass 2.
   !!
   !! **MO tensors here are physicist-ordered, `<pq|rs>`.** The gradient's AO
   !! side is chemist because the integrals are; these routines transcribe
   !! pycc's `_skeleton_lagrangian` and `_augment_with_canonical_pair_rotations`
   !! index for index, and holding the MO tensors in that code's order is what
   !! makes the transcription checkable line against line rather than through a
   !! layer of swaps. `<pq|rs> = (pr|qs)`; the reorder happens once, where the
   !! AO transform lands (`tools/mp2_hessian_oracle/CONVENTIONS.md` s.1-2).
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, atom_ao_blocks
   use mqc_libcint_hess_ints, only: hess_1e_block, hess_rinv_block, &
                                    hess_2e_skeleton_contract, eri_ip1_block, &
                                    HESS_KIN_II, HESS_KIN_IJ, HESS_NUC_II, HESS_NUC_IJ, &
                                    HESS_OVLP_II, HESS_OVLP_IJ, HESS_RINV_II, HESS_RINV_IJ
   use mqc_libcint_hessian, only: hcore_deriv_atom, overlap_deriv_atom
   implicit none
   private

   public :: mp2_skeleton_hessian
   public :: mp2_first_order_skeletons
   public :: mp2_mo_eri_physicist
   public :: mp2_cumulant_2pdm
   public :: mp2_mo_lagrangian
   public :: mp2_skeleton_lagrangian
   public :: mp2_pair_rotation_augment

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
      integer, allocatable :: owner(:), offsets(:), counts(:)
      real(dp) :: total_c(3, 3), total_r(3, 3)
      real(dp) :: wc, wr, zc
      integer :: n_ao, n_mo, natm, iao, jao, ia, ja
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
      ! Quartet-driven: evaluate, deposit, forget. The general orbit-weight
      ! fold lives with the integrals in `hess_2e_skeleton_contract`, and the
      ! reference deposits from the same buffers, so the second-derivative
      ! ERIs -- the dominant cost -- are generated exactly once.
      call hess_2e_skeleton_contract(mol, gamma_eff_ao, dref, owner, &
                                     hess_corr, hess_ref, error)
      if (error%has_error()) return

      deallocate (drel_ao, wcorr, wref, dref, owner, offsets, counts)
   end subroutine mp2_skeleton_hessian

   subroutine ao_to_mo_chem(a_ao, coeff, a_mo)
      !! Four-index MO transform of a chemist-ordered AO tensor
      !!
      !! Dense and unblocked on purpose: every Phase 1 tensor fits with room
      !! to spare (s.5 of the conventions note), and the blocked analog
      !! belongs to the unit that first needs it, not here.
      real(dp), intent(in) :: a_ao(:, :, :, :)
      real(dp), intent(in) :: coeff(:, :)
      real(dp), allocatable, intent(out) :: a_mo(:, :, :, :)

      real(dp), allocatable :: half(:, :, :, :)
      integer :: n_ao, n_mo, la, si, p, q

      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)
      allocate (half(n_mo, n_mo, n_ao, n_ao), a_mo(n_mo, n_mo, n_mo, n_mo))
      do si = 1, n_ao
         do la = 1, n_ao
            half(:, :, la, si) = matmul(transpose(coeff), &
                                        matmul(a_ao(:, :, la, si), coeff))
         end do
      end do
      do q = 1, n_mo
         do p = 1, n_mo
            a_mo(p, q, :, :) = matmul(transpose(coeff), &
                                      matmul(half(p, q, :, :), coeff))
         end do
      end do
      deallocate (half)
   end subroutine ao_to_mo_chem

   subroutine mp2_first_order_skeletons(mol, coeff, n_occ, fx, sx, erix, error)
      !! Unit 1.4: the per-perturbation first-order skeleton derivatives
      !!
      !! For every nuclear perturbation `X = (atom, component)`, at **fixed**
      !! MO coefficients:
      !!
      !!     S^(X)_pq   = C^T dS/dX C
      !!     f^(X)_pq   = C^T dh/dX C + sum_m [2 <pm|qm>^(X) - <pm|mq>^(X)]
      !!     <pq|rs>^(X)
      !!
      !! with `m` over the **full** occupied space, core included -- pycc's
      !! `fX` (`_hessian_blocks`), whose trace convention frozen-core work
      !! inherits. The perturbations stack as `x = 3*(atom-1) + component`,
      !! matching the oracle's `[3*atom + cart]`.
      !!
      !! The AO derivative is assembled the way `make_h1_atom` does it: the
      !! library puts the nabla on the bra's first index, so each of the four
      !! positions is permuted into first place in turn and kept when its
      !! function sits on the perturbed atom, with the library's sign.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)            !! C, (n_ao, n_mo)
      integer, intent(in) :: n_occ                   !! Full doubly-occupied count
      real(dp), allocatable, intent(out) :: fx(:, :, :)   !! (n_mo, n_mo, 3*natm)
      real(dp), allocatable, intent(out) :: sx(:, :, :)   !! (n_mo, n_mo, 3*natm)
      real(dp), allocatable, intent(out) :: erix(:, :, :, :, :)
         !! `<pq|rs>^(X)`, physicist, (n_mo, n_mo, n_mo, n_mo, 3*natm)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ip1(:, :, :, :, :), dao(:, :, :, :), chem(:, :, :, :)
      real(dp), allocatable :: h1(:, :, :), s1(:, :, :)
      integer, allocatable :: owner(:), offsets(:), counts(:)
      real(dp) :: acc
      integer :: n_ao, n_mo, natm, a, comp, x, c
      integer :: mu, nu, la, si, p, q, r, s, m

      if (error%has_error()) return

      n_ao = mol%nao
      n_mo = size(coeff, 2)
      natm = mol%natm
      if (size(coeff, 1) /= n_ao .or. n_occ < 1 .or. n_occ > n_mo) then
         call error%set(ERROR_VALIDATION, "the first-order skeletons were asked "// &
                        "for with coefficients that do not fit the basis.")
         return
      end if

      call eri_ip1_block(mol, ip1, error)
      if (error%has_error()) return

      allocate (offsets(natm), counts(natm), owner(n_ao))
      call atom_ao_blocks(mol, offsets, counts)
      owner = 0
      do c = 1, natm
         do mu = offsets(c) + 1, offsets(c) + counts(c)
            owner(mu) = c
         end do
      end do

      allocate (fx(n_mo, n_mo, 3*natm), sx(n_mo, n_mo, 3*natm))
      allocate (erix(n_mo, n_mo, n_mo, n_mo, 3*natm))
      allocate (dao(n_ao, n_ao, n_ao, n_ao))

      do a = 1, natm
         call hcore_deriv_atom(mol, a, h1, error)
         call overlap_deriv_atom(mol, a, s1, error)
         if (error%has_error()) return

         do comp = 1, 3
            x = 3*(a - 1) + comp
            sx(:, :, x) = matmul(transpose(coeff), matmul(s1(:, :, comp), coeff))
            fx(:, :, x) = matmul(transpose(coeff), matmul(h1(:, :, comp), coeff))

            dao = 0.0_dp
            do si = 1, n_ao
               do la = 1, n_ao
                  do nu = 1, n_ao
                     do mu = 1, n_ao
                        if (owner(mu) == a) dao(mu, nu, la, si) = &
                           dao(mu, nu, la, si) - ip1(mu, nu, la, si, comp)
                        if (owner(nu) == a) dao(mu, nu, la, si) = &
                           dao(mu, nu, la, si) - ip1(nu, mu, la, si, comp)
                        if (owner(la) == a) dao(mu, nu, la, si) = &
                           dao(mu, nu, la, si) - ip1(la, si, mu, nu, comp)
                        if (owner(si) == a) dao(mu, nu, la, si) = &
                           dao(mu, nu, la, si) - ip1(si, la, mu, nu, comp)
                     end do
                  end do
               end do
            end do

            call ao_to_mo_chem(dao, coeff, chem)

            ! The Fock skeleton's two-electron half, traced over the full
            ! occupied space: L^(X) folded on the spot rather than stored.
            do q = 1, n_mo
               do p = 1, n_mo
                  acc = 0.0_dp
                  do m = 1, n_occ
                     acc = acc + 2.0_dp*chem(p, q, m, m) - chem(p, m, m, q)
                  end do
                  fx(p, q, x) = fx(p, q, x) + acc
               end do
            end do

            ! Chemist to physicist, `<pq|rs> = (pr|qs)`, once and here only.
            do s = 1, n_mo
               do r = 1, n_mo
                  do q = 1, n_mo
                     do p = 1, n_mo
                        erix(p, q, r, s, x) = chem(p, r, q, s)
                     end do
                  end do
               end do
            end do
            deallocate (chem)
         end do
         deallocate (h1, s1)
      end do

      deallocate (ip1, dao, owner, offsets, counts)
   end subroutine mp2_first_order_skeletons

   subroutine mp2_mo_eri_physicist(mol, coeff, eri_mo, error)
      !! The unperturbed MO integrals, physicist order, `<pq|rs> = (pr|qs)`
      !!
      !! The unperturbed companion of `mp2_first_order_skeletons`' `erix`:
      !! what the Lagrangian and the pair-rotation augmentation contract.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)
      real(dp), allocatable, intent(out) :: eri_mo(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: ao(:, :, :, :), chem(:, :, :, :)
      integer :: n_mo, p, q, r, s

      if (error%has_error()) return
      if (size(coeff, 1) /= mol%nao) then
         call error%set(ERROR_VALIDATION, "the MO integral transform was handed "// &
                        "coefficients that do not fit the basis.")
         return
      end if

      call mol%eris(ao)
      call ao_to_mo_chem(ao, coeff, chem)
      deallocate (ao)

      n_mo = size(coeff, 2)
      allocate (eri_mo(n_mo, n_mo, n_mo, n_mo))
      do s = 1, n_mo
         do r = 1, n_mo
            do q = 1, n_mo
               do p = 1, n_mo
                  eri_mo(p, q, r, s) = chem(p, r, q, s)
               end do
            end do
         end do
      end do
      deallocate (chem)
   end subroutine mp2_mo_eri_physicist

   subroutine mp2_cumulant_2pdm(t2, n_frozen, n_occ, n_mo, gam)
      !! pycc's cumulant two-particle density, full MO, physicist order
      !!
      !!     Gamma_ijab = 2 t2_ijab - t2_ijba,   Gamma_abij = Gamma_ijab
      !!
      !! Only the oovv/vvoo blocks are nonzero for MP2 (`_mp2_tpdm`). This is
      !! **not** the gradient's `gamma_ao` and must not be conflated with it:
      !! that one carries the energy's coefficient and a different
      !! symmetrisation (conventions note, s.4a). Active amplitudes land at
      !! the active slices, so frozen rows and columns stay zero.
      real(dp), intent(in) :: t2(:, :, :, :)   !! (i, j, a, b), active occupied
      integer, intent(in) :: n_frozen, n_occ, n_mo
      real(dp), allocatable, intent(out) :: gam(:, :, :, :)

      real(dp) :: u
      integer :: i, j, a, b, n_o, n_v

      n_o = size(t2, 1)
      n_v = size(t2, 3)
      allocate (gam(n_mo, n_mo, n_mo, n_mo))
      gam = 0.0_dp
      do b = 1, n_v
         do a = 1, n_v
            do j = 1, n_o
               do i = 1, n_o
                  u = 2.0_dp*t2(i, j, a, b) - t2(i, j, b, a)
                  gam(n_frozen + i, n_frozen + j, n_occ + a, n_occ + b) = u
                  gam(n_occ + a, n_occ + b, n_frozen + i, n_frozen + j) = u
               end do
            end do
         end do
      end do
   end subroutine mp2_cumulant_2pdm

   subroutine mp2_mo_lagrangian(eri_mo, orbital_energies, drel, gam, n_occ, imat)
      !! The generalized-Fock orbital Lagrangian `I'`, pycc's `I`
      !!
      !!     I'_pq = -1/2 [ eps_p (D_pq + D_qp)
      !!                    + delta_{q occ} D_rs (L_rpsq + L_rqsp)
      !!                    + 4 <pr|st> Gamma_qrst ]
      !!
      !! Transcribed from `_spatial_lagrangian`. This is the correlation-only
      !! matrix that multiplies `dS/dR`; the gradient's `imat` is a *part* of
      !! it and its `energy_weighted_ao` is the total including the
      !! reference's share -- the reconciliations are in the conventions note,
      !! s.4a, and `test_mqc_mp2_hessian_response` checks this routine against
      !! that reconstruction.
      real(dp), intent(in) :: eri_mo(:, :, :, :)   !! Physicist, unperturbed
      real(dp), intent(in) :: orbital_energies(:)
      real(dp), intent(in) :: drel(:, :)           !! Relaxed 1-PDM, MO
      real(dp), intent(in) :: gam(:, :, :, :)      !! From `mp2_cumulant_2pdm`
      integer, intent(in) :: n_occ                 !! Full occupied count
      real(dp), allocatable, intent(out) :: imat(:, :)

      real(dp) :: acc
      integer :: n_mo, p, q, r, s, t

      n_mo = size(drel, 1)
      allocate (imat(n_mo, n_mo))

      do q = 1, n_mo
         do p = 1, n_mo
            acc = orbital_energies(p)*(drel(p, q) + drel(q, p))
            do t = 1, n_mo
               do s = 1, n_mo
                  do r = 1, n_mo
                     acc = acc + 4.0_dp*eri_mo(p, r, s, t)*gam(q, r, s, t)
                  end do
               end do
            end do
            imat(p, q) = acc
         end do
      end do
      do q = 1, n_occ
         do p = 1, n_mo
            acc = 0.0_dp
            do s = 1, n_mo
               do r = 1, n_mo
                  acc = acc + drel(r, s)* &
                        (2.0_dp*eri_mo(r, p, s, q) - eri_mo(r, p, q, s) &
                         + 2.0_dp*eri_mo(r, q, s, p) - eri_mo(r, q, p, s))
               end do
            end do
            imat(p, q) = imat(p, q) + acc
         end do
      end do
      imat = -0.5_dp*imat
   end subroutine mp2_mo_lagrangian

   subroutine mp2_skeleton_lagrangian(fx1, sx1, erix1, drel, gam, imat, n_occ, &
                                      ip, xov, i2)
      !! Unit 1.5: the skeleton-perturbed orbital Lagrangian for one `X`
      !!
      !! `_skeleton_lagrangian`, transcribed at fixed unperturbed densities:
      !!
      !!     I'^(X)_pq = -1/2 [ D_qr f^(X)_pr + D_rq f^(X)_rp
      !!                        + delta_{q occ} D_rs (L^(X)_rpsq + L^(X)_rqsp)
      !!                        + 4 <pr|st>^(X) Gamma_qrst
      !!                        + I_qr S^(X)_pr + I_rq S^(X)_rp ]
      !!
      !! with `X^(X) = I'^(X)T - I'^(X)` the orbital-response driver and
      !! `I''^(X)` the energy-weighted rewrite -- `I'^(X)` with its
      !! virtual-occupied block replaced by the transposed occupied-virtual
      !! one. All-electron non-canonical work uses these as they are; a frozen
      !! core or the canonical gauge sends them through
      !! `mp2_pair_rotation_augment` first.
      real(dp), intent(in) :: fx1(:, :), sx1(:, :)   !! One perturbation's f^(X), S^(X)
      real(dp), intent(in) :: erix1(:, :, :, :)      !! `<pq|rs>^(X)`, physicist
      real(dp), intent(in) :: drel(:, :)             !! Relaxed 1-PDM, MO
      real(dp), intent(in) :: gam(:, :, :, :)        !! From `mp2_cumulant_2pdm`
      real(dp), intent(in) :: imat(:, :)             !! From `mp2_mo_lagrangian`
      integer, intent(in) :: n_occ                   !! Full occupied count
      real(dp), allocatable, intent(out) :: ip(:, :)
      real(dp), allocatable, intent(out) :: xov(:, :)
      real(dp), allocatable, intent(out) :: i2(:, :)

      real(dp) :: acc
      integer :: n_mo, p, q, r, s, t

      n_mo = size(drel, 1)
      allocate (ip(n_mo, n_mo), xov(n_mo, n_mo), i2(n_mo, n_mo))

      ! One- and two-particle terms with both derivative placements; the
      ! occupied-column two-electron term folds L^(X) on the spot, as the
      ! Fock skeleton build does.
      do q = 1, n_mo
         do p = 1, n_mo
            acc = 0.0_dp
            do r = 1, n_mo
               acc = acc + drel(q, r)*fx1(p, r) + drel(r, q)*fx1(r, p) &
                     + imat(q, r)*sx1(p, r) + imat(r, q)*sx1(r, p)
            end do
            do t = 1, n_mo
               do s = 1, n_mo
                  do r = 1, n_mo
                     acc = acc + 4.0_dp*erix1(p, r, s, t)*gam(q, r, s, t)
                  end do
               end do
            end do
            ip(p, q) = acc
         end do
      end do
      do q = 1, n_occ
         do p = 1, n_mo
            acc = 0.0_dp
            do s = 1, n_mo
               do r = 1, n_mo
                  acc = acc + drel(r, s)* &
                        (2.0_dp*erix1(r, p, s, q) - erix1(r, p, q, s) &
                         + 2.0_dp*erix1(r, q, s, p) - erix1(r, q, p, s))
               end do
            end do
            ip(p, q) = ip(p, q) + acc
         end do
      end do
      ip = -0.5_dp*ip

      xov = transpose(ip) - ip
      i2 = ip
      i2(n_occ + 1:n_mo, 1:n_occ) = transpose(ip(1:n_occ, n_occ + 1:n_mo))
   end subroutine mp2_skeleton_lagrangian

   subroutine mp2_pair_rotation_augment(ip, xov, i2, l_mo, orbital_energies, &
                                        n_occ, n_frozen, canonical, xt, it, pf)
      !! The closed-form pair rotations folded into the skeleton carriers
      !!
      !! `_augment_with_canonical_pair_rotations`, transcribed. Two rotations
      !! the occupied-virtual solve does not provide, both from the divide
      !! `P^(X)_pq = (I'^(X)_pq - I'^(X)_qp) / (eps_p - eps_q)`:
      !!
      !! * the **independent** core<->active-occupied block, always present
      !!   with a frozen core -- an ungated direct divide, its gap never small;
      !! * the **redundant** active oo/vv blocks, canonical gauge only --
      !!   gap-gated at 1e-8 (`_dependent_pairs`), because near a degeneracy
      !!   the numerator vanishes by the same symmetry that makes the divide
      !!   ill-conditioned.
      !!
      !! All-electron non-canonical, this is arithmetically the identity --
      !! both guards skip and `xt`, `it` copy through. It exists now, written
      !! against pycc while the structure is fresh, because Phase 2's frozen
      !! core reaches it (the plan's Unit 1.5 note); nothing in Phase 1 does.
      real(dp), intent(in) :: ip(:, :), xov(:, :), i2(:, :)
      real(dp), intent(in) :: l_mo(:, :, :, :)
         !! Unperturbed orbital-Hessian weight `L_pqrs = 2 <pq|rs> - <pq|sr>`
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ                   !! Full occupied count
      integer, intent(in) :: n_frozen
      logical, intent(in) :: canonical
      real(dp), allocatable, intent(out) :: xt(:, :)
      real(dp), allocatable, intent(out) :: it(:, :)
      real(dp), allocatable, intent(out) :: pf(:, :)

      real(dp), allocatable :: pof(:, :), pvv(:, :), aug(:, :)
      real(dp) :: acc
      integer :: n_mo, n_vir, p, q, k, l, a, i

      n_mo = size(ip, 1)
      n_vir = n_mo - n_occ
      allocate (pof(n_occ, n_occ), pvv(n_vir, n_vir), pf(n_mo, n_mo))
      pof = 0.0_dp
      pvv = 0.0_dp

      if (n_frozen > 0) then
         ! Independent core<->active-occupied rotation: rows core, columns
         ! active, then mirrored -- `P` is symmetric.
         do q = n_frozen + 1, n_occ
            do p = 1, n_frozen
               pof(p, q) = (ip(p, q) - ip(q, p)) &
                           /(orbital_energies(p) - orbital_energies(q))
               pof(q, p) = pof(p, q)
            end do
         end do
      end if
      if (canonical) then
         call dependent_pairs(ip(n_frozen + 1:n_occ, n_frozen + 1:n_occ), &
                              orbital_energies(n_frozen + 1:n_occ), aug)
         pof(n_frozen + 1:n_occ, n_frozen + 1:n_occ) = aug
         deallocate (aug)
         call dependent_pairs(ip(n_occ + 1:n_mo, n_occ + 1:n_mo), &
                              orbital_energies(n_occ + 1:n_mo), aug)
         pvv = aug
         deallocate (aug)
      end if

      pf = 0.0_dp
      pf(1:n_occ, 1:n_occ) = pof
      pf(n_occ + 1:n_mo, n_occ + 1:n_mo) = pvv

      ! X~^(X)_ai: the occupied pair sums run over the full occupied space.
      xt = xov
      do i = 1, n_occ
         do a = n_occ + 1, n_mo
            acc = 0.0_dp
            do l = 1, n_occ
               do k = 1, n_occ
                  acc = acc + pof(k, l)*(l_mo(k, a, l, i) + l_mo(k, i, l, a))
               end do
            end do
            do l = 1, n_vir
               do k = 1, n_vir
                  acc = acc + pvv(k, l)*(l_mo(n_occ + k, a, n_occ + l, i) &
                                         + l_mo(n_occ + k, i, n_occ + l, a))
               end do
            end do
            xt(a, i) = xt(a, i) + 0.5_dp*acc
         end do
      end do

      ! I~''^(X): the occupied block picks up the rotation's energy weight and
      ! the same pair sums; the virtual block only the energy weight.
      it = i2
      do q = 1, n_occ
         do p = 1, n_occ
            acc = 0.0_dp
            do l = 1, n_occ
               do k = 1, n_occ
                  acc = acc + pof(k, l)*(l_mo(k, p, l, q) + l_mo(k, q, l, p))
               end do
            end do
            do l = 1, n_vir
               do k = 1, n_vir
                  acc = acc + pvv(k, l)*(l_mo(n_occ + k, p, n_occ + l, q) &
                                         + l_mo(n_occ + k, q, n_occ + l, p))
               end do
            end do
            it(p, q) = it(p, q) - pof(p, q)*orbital_energies(q) - 0.5_dp*acc
         end do
      end do
      do q = 1, n_vir
         do p = 1, n_vir
            it(n_occ + p, n_occ + q) = it(n_occ + p, n_occ + q) &
                                       - pvv(p, q)*orbital_energies(n_occ + q)
         end do
      end do

      deallocate (pof, pvv)
   end subroutine mp2_pair_rotation_augment

   subroutine dependent_pairs(iblock, eps_block, p)
      !! `P_mn = (I_mn - I_nm) / (eps_m - eps_n)`, gap-gated
      !!
      !! The gate is on the **gap**, not the numerator: a numerator gate would
      !! zero small-but-nonzero rotations whose derivative is not zero,
      !! leaving `P` and `dP` inconsistent -- pycc's `_dependent_pairs`
      !! records that as a genuine low-symmetry Hessian error, not a style
      !! choice.
      real(dp), intent(in) :: iblock(:, :)
      real(dp), intent(in) :: eps_block(:)
      real(dp), allocatable, intent(out) :: p(:, :)

      real(dp), parameter :: GAP_THRESH = 1.0e-8_dp
      real(dp) :: gap
      integer :: m, n

      allocate (p(size(iblock, 1), size(iblock, 2)))
      p = 0.0_dp
      do n = 1, size(iblock, 2)
         do m = 1, size(iblock, 1)
            gap = eps_block(m) - eps_block(n)
            if (abs(gap) > GAP_THRESH) then
               p(m, n) = (iblock(m, n) - iblock(n, m))/gap
            end if
         end do
      end do
   end subroutine dependent_pairs

end module mqc_libcint_mp2_hessian
