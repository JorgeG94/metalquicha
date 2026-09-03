!! The fixed-density second-derivative skeleton of the MP2 correlation Hessian
module mqc_czt_mp2_hessian
   !! The analytic MP2 correlation Hessian, in two passes, transcribed from
   !! pycc's `_hessian_blocks` (`route == 'aod'`).
   !!
   !! Pass 1, `mp2_skeleton_hessian`: the per-atom-pair scalar
   !!
   !!     s = Gamma_eff . (mu nu|la si)^{XY} + D_rel . h^{XY} + I . S^{XY}
   !!
   !! contracted against the **unperturbed** densities. The SCF reference's
   !! skeleton rides the same sweep -- the second-derivative two-electron
   !! integrals dominate the cost and both blocks contract the same ones -- so
   !! they are generated once and deposited into two accumulators. That also
   !! makes the reference block checkable against `partial_hessian`, which
   !! walks its own quartets.
   !!
   !! Pass 2: the per-perturbation first-order skeletons `f^(X)`, `S^(X)`,
   !! `<pq|rs>^(X)`; the skeleton Lagrangian `I'^(X)` and its orbital-response
   !! carriers; and the perturbed response, where the CPHF-folded `d_Y f` and
   !! `d_Y <pq|rs>` feed the closed-form perturbed amplitudes, and those feed
   !! the perturbed relaxed density, energy-weighted density and cumulant
   !! through one Z-vector solve per perturbation.
   !!
   !! **Storage.** The per-perturbation `d_Y <pq|rs>` and its spin-adapted `L`
   !! live one at a time inside the response driver's loop; what the perturbed
   !! energy-weighted density needs after the batched Z-vector solve is folded
   !! out *before* the tensor is dropped, so nothing `nmo^4` is ever
   !! `3N`-resident. The exception is the skeleton stack `erix`, which is dense
   !! `nmo^4 x 3N`.
   !!
   !! **MO tensors here are physicist-ordered, `<pq|rs>`**, matching the pycc
   !! routines these transcribe index for index, where the gradient's AO side
   !! is chemist. `<pq|rs> = (pr|qs)`; the reorder happens once, where the AO
   !! transform lands (`tools/mp2_hessian_oracle/CONVENTIONS.md` s.1-2).
   !!
   !! **The same routines serve the double hybrid.** With `xc` present, every
   !! reference-operator weight becomes `2 <pq|rs> - k <pq|sr>` plus an
   !! exchange-correlation kernel term. No four-index `f_xc` tensor is built for
   !! it: every kernel contraction is the reference operator applied to a
   !! symmetric generalized density, through `ks_pair_kernel` and its derivative
   !! companions at the bottom of this file.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks
   use mqc_czt_hess_ints, only: hess_1e_block, hess_rinv_block, &
                                hess_2e_skeleton_contract, eri_ip1_block, &
                                HESS_KIN_II, HESS_KIN_IJ, HESS_NUC_II, HESS_NUC_IJ, &
                                HESS_OVLP_II, HESS_OVLP_IJ, HESS_RINV_II, HESS_RINV_IJ
   use mqc_czt_hessian, only: hcore_deriv_atom, overlap_deriv_atom, &
                              make_h1_atom, solve_mo1_batch
   use mqc_czt_cphf, only: cphf_solve
   use mqc_czt_mp2, only: transform_ovov
   use mqc_czt_mp2_gradient, only: czt_mp2_gradient, build_amplitudes, &
                                   build_effective_2pdm_ao
   use mqc_czt_multipole, only: multipole_matrices, dipole_integral_derivatives
   use mqc_czt_xc, only: xc_context_t, xc_kernel_apply, xc_kernel2_apply
   use mqc_czt_xc_hessian, only: xc_potential_deriv, xc_kernel_deriv, &
                                 xc_potential_hessian
   implicit none
   private

   public :: mp2_skeleton_hessian
   public :: mp2_first_order_skeletons
   public :: mp2_mo_eri_physicist
   public :: mp2_cumulant_2pdm
   public :: mp2_mo_lagrangian
   public :: mp2_skeleton_lagrangian
   public :: mp2_pair_rotation_augment
   public :: mp2_orbital_response_term
   public :: mp2_full_u
   public :: mp2_perturbed_fock
   public :: mp2_perturbed_eri
   public :: mp2_perturbed_t2
   public :: mp2_perturbed_onepdm
   public :: mp2_perturbed_lagrangian
   public :: mp2_perturbed_zvector_rhs
   public :: mp2_perturbed_response
   public :: mp2_correlation_hessian
   public :: mp2_ks_first_order_fock

contains

   subroutine mp2_skeleton_hessian(mol, gamma_eff_ao, relaxed_density_mo, &
                                   energy_weighted_ao, coeff, orbital_energies, &
                                   n_occ, hess_corr, hess_ref, error, xc, &
                                   scf_density)
      !! Both fixed-density skeletons, `(3, 3, natm, natm)` each
      !!
      !! `hess_corr(a, b, A, B)` is pycc's pass-1 scalar `s` for the atom pair
      !! `(A, B)` and Cartesian components `(a, b)`; `hess_ref` is the SCF
      !! skeleton over the same integrals -- the `partial_hessian` numbers,
      !! deposited from this sweep rather than a second one.
      !!
      !! **`energy_weighted_ao` is the total, not the correlation share.** The
      !! gradient's `W` carries the reference's energy-weighted density, which
      !! is most of its norm, where pycc's `I` is correlation only. The split
      !! is done here, once:
      !!
      !!     I_corr = (W_total + 2 W_ref) / 2,   W_ref = 2 sum_i eps_i C_i C_i^T
      !!
      !! (`tools/mp2_hessian_oracle/CONVENTIONS.md` s.4a). Handing this routine
      !! a correlation-only `W` double-counts the reference's share -- give it
      !! exactly what `czt_mp2_gradient` returns.
      !!
      !! **With `xc` the convention flips.** The double-hybrid gradient's `W`
      !! never took the reference's share on board, so here it is halved and
      !! nothing is added back. The sweep changes twice more: the reference
      !! deposit's exchange carries the functional's fraction, and
      !! `Tr(D_rel V_xc^{XY})` -- the piece of `D_rel f^{XY}` no density fold
      !! can carry -- lands from `xc_potential_hessian`. `hess_ref` is still the
      !! Hartree-Fock-shaped skeleton over `dref` and is to be **discarded** by
      !! a double-hybrid caller, whose reference block comes whole from
      !! `ks_hessian`.
      type(czt_molecule_t), intent(in) :: mol
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
      type(xc_context_t), intent(inout), optional :: xc
         !! The functional the reference used. Present switches this sweep to
         !! the double hybrid's conventions -- see the header.
      real(dp), intent(in), optional :: scf_density(:, :)
         !! The converged reference density, where `V_xc` is evaluated.
         !! Required with `xc`.

      real(dp), allocatable :: drel_ao(:, :), wcorr(:, :), wref(:, :), dref(:, :)
      real(dp), allocatable :: h2(:, :, :), hab(:, :, :), s2(:, :, :), sab(:, :, :)
      real(dp), allocatable :: tmp(:, :, :), r2(:, :, :), rab(:, :, :)
      real(dp), allocatable :: cross_c(:, :, :), cross_r(:, :, :)
      integer, allocatable :: owner(:), offsets(:), counts(:)
      real(dp) :: total_c(3, 3), total_r(3, 3)
      real(dp) :: wc, wr, zc, kf
      logical :: dh
      integer :: n_ao, n_mo, natm, iao, jao, ia, ja
      integer :: a, b, comp, c, i

      if (error%has_error()) return

      dh = present(xc)
      kf = 1.0_dp
      if (dh) then
         if (.not. present(scf_density)) then
            call error%set(ERROR_VALIDATION, "the double-hybrid skeleton Hessian "// &
                           "needs the converged reference density alongside the "// &
                           "functional context")
            return
         end if
         kf = xc%exx_fraction
      end if

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
      if (dh) then
         ! The double-hybrid gradient's `W` never carried the reference share,
         ! so there is nothing to strip back out (header).
         wcorr = 0.5_dp*energy_weighted_ao
      else
         wcorr = 0.5_dp*(energy_weighted_ao + 2.0_dp*wref)
      end if

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
      ! The factor of two is bra and ket contributing the same number; the signs
      ! follow the energy expression -- `+ D h^{XY}` against `- W_ref S^{XY}`
      ! for the reference, while the correlation's energy-weighted term enters
      ! with pycc's sign, `+ I S^{XY}`.
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
      ! term for term. Second-derivative one-electron integrals on a steep core
      ! exponent reach ~1e5, so the accumulator swings that far before
      ! cancelling to order one, and any regrouping -- one fused update instead
      ! of an add and a subtract -- moves the reference block by ~1e-11, which
      ! is what the 1e-13 gate against `partial_hessian` catches. Do not fuse
      ! them back for tidiness.
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
      ! The operator-centre terms of the nuclear attraction, as
      ! `partial_hessian` assembles them: once per nucleus, both densities
      ! deposited from the same integral blocks.
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
      ! Quartet-driven: evaluate, deposit, forget. The orbit-weight fold lives
      ! with the integrals in `hess_2e_skeleton_contract`, and the reference
      ! deposits from the same buffers, so the second-derivative ERIs -- the
      ! dominant cost -- are generated exactly once.
      call hess_2e_skeleton_contract(mol, gamma_eff_ao, dref, owner, &
                                     hess_corr, hess_ref, error, k_scale=kf)
      if (error%has_error()) return

      ! The reference operator's second derivative contains the
      ! exchange-correlation potential's, contracted against the relaxed
      ! correlation density -- the one pass-1 term with no Hartree-Fock
      ! counterpart, since `gamma_eff_ao`'s fold can only carry the separable
      ! two-electron half of `D_rel f^{XY}`. Grid held fixed, as everywhere
      ! on this Hessian's ladder.
      if (dh) then
         call xc_potential_hessian(xc, mol, scf_density, drel_ao, hess_corr, error)
         if (error%has_error()) return
      end if

      deallocate (drel_ao, wcorr, wref, dref, owner, offsets, counts)
   end subroutine mp2_skeleton_hessian

   subroutine ao_to_mo_chem(a_ao, coeff, a_mo)
      !! Four-index MO transform of a chemist-ordered AO tensor
      !!
      !! Dense and unblocked -- every tensor on this path fits (s.5 of the
      !! conventions note).
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
      !! The per-perturbation first-order skeleton derivatives
      !!
      !! For every nuclear perturbation `X = (atom, component)`, at **fixed**
      !! MO coefficients:
      !!
      !!     S^(X)_pq   = C^T dS/dX C
      !!     f^(X)_pq   = C^T dh/dX C + sum_m [2 <pm|qm>^(X) - <pm|mq>^(X)]
      !!     <pq|rs>^(X)
      !!
      !! with `m` over the **full** occupied space, core included -- pycc's `fX`
      !! trace convention, which frozen-core work inherits. The perturbations
      !! stack as `x = 3*(atom-1) + component`, matching the oracle's
      !! `[3*atom + cart]`.
      !!
      !! The AO derivative is assembled the way `make_h1_atom` does it: the
      !! library puts the nabla on the bra's first index, so each of the four
      !! positions is permuted into first place in turn and kept when its
      !! function sits on the perturbed atom, with the library's sign.
      type(czt_molecule_t), intent(in) :: mol
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
      !! The unperturbed companion of `mp2_first_order_skeletons`' `erix`, and
      !! what the Lagrangian and the pair-rotation augmentation contract.
      type(czt_molecule_t), intent(in) :: mol
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
      !! Only the oovv/vvoo blocks are nonzero for MP2. This is **not** the
      !! gradient's `gamma_ao`: that one carries the energy's coefficient and a
      !! different symmetrisation (conventions note, s.4a). Active amplitudes
      !! land at the active slices, so frozen rows and columns stay zero.
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
      !! Transcribed from `_spatial_lagrangian`. The correlation-only matrix
      !! that multiplies `dS/dR`; the gradient's `imat` is a *part* of it and
      !! its `energy_weighted_ao` is the total including the reference's share
      !! (conventions note, s.4a).
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
      !! The skeleton-perturbed orbital Lagrangian for one `X`
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
      !! virtual-occupied block replaced by the transposed occupied-virtual one.
      !! All-electron non-canonical work uses these as they are; a frozen core
      !! or the canonical gauge sends them through `mp2_pair_rotation_augment`
      !! first.
      real(dp), intent(in) :: fx1(:, :), sx1(:, :)   !! One perturbation's f^(X), S^(X)
      real(dp), intent(in), target, contiguous :: erix1(:, :, :, :)
         !! `<pq|rs>^(X)`, physicist
      real(dp), intent(in) :: drel(:, :)             !! Relaxed 1-PDM, MO
      real(dp), intent(in), target, contiguous :: gam(:, :, :, :)
         !! From `mp2_cumulant_2pdm`
      real(dp), intent(in) :: imat(:, :)             !! From `mp2_mo_lagrangian`
      integer, intent(in) :: n_occ                   !! Full occupied count
      real(dp), allocatable, intent(out) :: ip(:, :)
      real(dp), allocatable, intent(out) :: xov(:, :)
      real(dp), allocatable, intent(out) :: i2(:, :)

      real(dp), pointer :: erix2(:, :), gam2(:, :)
      real(dp), allocatable :: two_e(:, :)
      real(dp) :: acc
      integer :: n_mo, n_cube, p, q, r, s

      n_mo = size(drel, 1)
      n_cube = n_mo*n_mo*n_mo
      allocate (ip(n_mo, n_mo), xov(n_mo, n_mo), i2(n_mo, n_mo), two_e(n_mo, n_mo))

      ! One- and two-particle terms with both derivative placements, each a
      ! matrix product. The four one-electron terms differ only in which operand
      ! is transposed, and the `4 <pr|st>^(X) Gamma_qrst` term contracts a
      ! contiguous trailing `(r,s,t)` against `p` and `q` leading, so it is one
      ! gemm rather than an `n_mo^5` nest.
      call pic_gemm(fx1, drel, ip, transb="T")
      call pic_gemm(fx1, drel, ip, transa="T", beta=1.0_dp)
      call pic_gemm(sx1, imat, ip, transb="T", beta=1.0_dp)
      call pic_gemm(sx1, imat, ip, transa="T", beta=1.0_dp)

      erix2(1:n_mo, 1:n_cube) => erix1
      gam2(1:n_mo, 1:n_cube) => gam
      call pic_gemm(erix2, gam2, two_e, transb="T")
      ip = ip + 4.0_dp*two_e

      ! The occupied-column term folds L^(X) on the spot, as the Fock skeleton
      ! build does. It stays a loop: `p` and `q` are strided inside `erix1`, and
      ! it is one power of `n_mo` cheaper than the term above, so the
      ! independent `(p,q)` is what gets parallelised.
      !$omp parallel do collapse(2) default(none) schedule(static) &
      !$omp    shared(ip, drel, erix1, n_mo, n_occ) private(p, q, r, s, acc)
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
      !$omp end parallel do
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
      !! All-electron non-canonical, this is arithmetically the identity: both
      !! guards skip and `xt`, `it` copy through.
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
      !! The gate is on the **gap**, not the numerator: a numerator gate zeroes
      !! small-but-nonzero rotations whose derivative is not zero, leaving `P`
      !! and `dP` inconsistent, which pycc's `_dependent_pairs` records as a
      !! genuine low-symmetry Hessian error.
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

   subroutine mp2_orbital_response_term(mo1, sx, fx, xt, i2t, orb, pf)
      !! The orbital-response share of the pass-2 response
      !!
      !!     orb(X, Y) = 2 U^Y_ai X~^(X)_ai + S^(Y)_pq I~''^(X)_pq
      !!                 [+ P^(X)_pq f^(Y)_pq when augmented]
      !!
      !! `U^Y` is `solve_mo1_batch`'s first-order orbitals -- virtual rows the
      !! solved response, occupied rows the `-S/2` orthonormality fixes -- and
      !! only its virtual-occupied block is read here. `pf` accompanies the
      !! augmented carriers and is meaningless without them.
      real(dp), intent(in) :: mo1(:, :, :, :)   !! (n_mo, n_occ, 3, natm)
      real(dp), intent(in) :: sx(:, :, :)       !! (n_mo, n_mo, 3*natm)
      real(dp), intent(in) :: fx(:, :, :)       !! Same layout; read only with `pf`
      real(dp), intent(in) :: xt(:, :, :)       !! X~^(X) (or bare X^(X))
      real(dp), intent(in) :: i2t(:, :, :)      !! I~''^(X) (or bare I''^(X))
      real(dp), allocatable, intent(out) :: orb(:, :)   !! (3*natm, 3*natm)
      real(dp), intent(in), optional :: pf(:, :, :)     !! P^(X), full MO

      real(dp) :: acc
      integer :: n_mo, n_occ, n_pert, ix, iy, ay, cy, a, i

      n_mo = size(mo1, 1)
      n_occ = size(mo1, 2)
      n_pert = size(sx, 3)
      allocate (orb(n_pert, n_pert))

      do iy = 1, n_pert
         ay = (iy - 1)/3 + 1
         cy = iy - 3*(ay - 1)
         do ix = 1, n_pert
            acc = 0.0_dp
            do i = 1, n_occ
               do a = n_occ + 1, n_mo
                  acc = acc + 2.0_dp*mo1(a, i, cy, ay)*xt(a, i, ix)
               end do
            end do
            acc = acc + sum(sx(:, :, iy)*i2t(:, :, ix))
            if (present(pf)) acc = acc + sum(pf(:, :, ix)*fx(:, :, iy))
            orb(ix, iy) = acc
         end do
      end do
   end subroutine mp2_orbital_response_term

   subroutine mp2_full_u(mo1_x, sx1, n_occ, u, n_frozen, fx1, l_mo, &
                         orbital_energies)
      !! The full `nmo x nmo` orbital rotation `U^X` for one perturbation
      !!
      !! pycc's `full_U`. `solve_mo1_batch`'s first-order orbitals already carry
      !! the whole occupied-column block -- the solved virtual rows and the
      !! `-S^(X)/2` occupied rows the orthonormality fixes -- so only the
      !! virtual columns are assembled here:
      !!
      !!     U^X_ab = -1/2 S^(X)_ab
      !!     U^X_ia = -S^(X)_ia - U^X_ai      (U + U^T = -S^(X))
      !!
      !! With a frozen core the core<->active-occupied rotation is
      !! **independent** -- it moves the frozen/active partition -- so that
      !! block is rewritten from the canonical Brillouin condition
      !! `d_X f_ci = 0` rather than left at the orthonormality value, with
      !! `U_ic` eliminated through `U_ic = -S^(X)_ci - U_ci`:
      !!
      !!     U^X_ci = -[ f^(X)_ci - S^(X)_ci eps_i
      !!                 - 1/2 S^(X)_nm (L_cnim + L_cmin)
      !!                 + U^X_dm (L_cdim + L_cmid) ] / (eps_c - eps_i)
      !!     U^X_ic = -S^(X)_ci - U^X_ci
      !!
      !! `n, m` over the full occupied space, `d` over the virtuals -- the
      !! virtual sum reads the solved response already sitting in `u`. The
      !! divide is ungated: the core-active gap is never small. The redundant
      !! core-core, active-active and virtual-virtual blocks stay at
      !! `-1/2 S^(X)`. `fx1`, `l_mo` and `orbital_energies` are read only
      !! when `n_frozen > 0`.
      real(dp), intent(in) :: mo1_x(:, :)   !! (n_mo, n_occ), one perturbation
      real(dp), intent(in) :: sx1(:, :)     !! S^(X), MO, (n_mo, n_mo)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: u(:, :)
      integer, intent(in), optional :: n_frozen
      real(dp), intent(in), optional :: fx1(:, :)   !! f^(X), MO, (n_mo, n_mo)
      real(dp), intent(in), optional :: l_mo(:, :, :, :)
         !! Unperturbed `L_pqrs = 2 <pq|rs> - <pq|sr>`
      real(dp), intent(in), optional :: orbital_energies(:)

      real(dp) :: acc
      integer :: n_mo, n_core, i, a, c, d, m, n

      n_mo = size(mo1_x, 1)
      n_core = 0
      if (present(n_frozen)) n_core = n_frozen
      allocate (u(n_mo, n_mo))
      u(:, 1:n_occ) = mo1_x
      u(n_occ + 1:n_mo, n_occ + 1:n_mo) = -0.5_dp*sx1(n_occ + 1:n_mo, n_occ + 1:n_mo)
      do a = n_occ + 1, n_mo
         do i = 1, n_occ
            u(i, a) = -sx1(i, a) - mo1_x(a, i)
         end do
      end do
      if (n_core > 0) then
         ! The Brillouin rewrite reads the solved virtual response out of the
         ! occupied columns just assembled, so it runs after them and touches
         ! only the core<->active sub-blocks.
         do i = n_core + 1, n_occ
            do c = 1, n_core
               acc = fx1(c, i) - sx1(c, i)*orbital_energies(i)
               do m = 1, n_occ
                  do n = 1, n_occ
                     acc = acc - 0.5_dp*sx1(n, m)*(l_mo(c, n, i, m) + l_mo(c, m, i, n))
                  end do
               end do
               do m = 1, n_occ
                  do d = n_occ + 1, n_mo
                     acc = acc + u(d, m)*(l_mo(c, d, i, m) + l_mo(c, m, i, d))
                  end do
               end do
               u(c, i) = -acc/(orbital_energies(c) - orbital_energies(i))
               u(i, c) = -sx1(c, i) - u(c, i)
            end do
         end do
      end if
   end subroutine mp2_full_u

   subroutine mp2_perturbed_fock(fx1, sx1, u, l_mo, orbital_energies, n_occ, df)
      !! The full first derivative of the MO Fock matrix, `d_X f`
      !!
      !! pycc's `CPHF.perturbed_fock`: the skeleton plus the orbital response,
      !!
      !!     d_X f_pq = f^(X)_pq + U^X_pq eps_p + U^X_qp eps_q
      !!                - 1/2 S^(X)_nm (L_pnqm + L_pmqn)
      !!                + U^X_cm (L_pcqm + L_pmqc)
      !!
      !! with `n, m` over the full occupied space and `c` over the virtuals. The
      !! oo/vv blocks are the non-canonical coupling the perturbed amplitudes
      !! fold in and the diagonal is `d_X eps`, but it is the **whole matrix**
      !! the perturbed Lagrangian multiplies, never just that diagonal.
      real(dp), intent(in) :: fx1(:, :), sx1(:, :)   !! f^(X), S^(X), MO
      real(dp), intent(in) :: u(:, :)                !! From `mp2_full_u`
      real(dp), intent(in) :: l_mo(:, :, :, :)
         !! Unperturbed `L_pqrs = 2 <pq|rs> - <pq|sr>`
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: df(:, :)

      real(dp) :: acc
      integer :: n_mo, p, q, m, n, c

      n_mo = size(fx1, 1)
      allocate (df(n_mo, n_mo))
      do q = 1, n_mo
         do p = 1, n_mo
            acc = fx1(p, q) + u(p, q)*orbital_energies(p) + u(q, p)*orbital_energies(q)
            do m = 1, n_occ
               do n = 1, n_occ
                  acc = acc - 0.5_dp*sx1(n, m)*(l_mo(p, n, q, m) + l_mo(p, m, q, n))
               end do
            end do
            do m = 1, n_occ
               do c = n_occ + 1, n_mo
                  acc = acc + u(c, m)*(l_mo(p, c, q, m) + l_mo(p, m, q, c))
               end do
            end do
            df(p, q) = acc
         end do
      end do
   end subroutine mp2_perturbed_fock

   subroutine mp2_perturbed_eri(erix1, u, eri_mo, deri)
      !! The full first derivative of the MO integrals, `d_X <pq|rs>`
      !!
      !! pycc's `CPHF.perturbed_eri`: the skeleton plus all four rotations,
      !!
      !!     d_X <pq|rs> = <pq|rs>^(X) + U^X_tp <tq|rs> + U^X_tq <pt|rs>
      !!                              + U^X_tr <pq|ts> + U^X_ts <pq|rt>
      !!
      !! Each rotation is one matrix product over a reshaped view; the caller
      !! holds one perturbation's tensor at a time (module header, storage).
      real(dp), intent(in) :: erix1(:, :, :, :)   !! `<pq|rs>^(X)`, physicist
      real(dp), intent(in) :: u(:, :)
      real(dp), intent(in) :: eri_mo(:, :, :, :)  !! Unperturbed, physicist
      real(dp), allocatable, intent(out) :: deri(:, :, :, :)

      integer :: n, r, s

      n = size(eri_mo, 1)
      allocate (deri(n, n, n, n))
      ! First index: (t, [qrs]) -> (p, [qrs]).
      deri = reshape(matmul(transpose(u), reshape(eri_mo, [n, n*n*n])), [n, n, n, n])
      ! Second index, per ket pair; third index, per ket component.
      do s = 1, n
         do r = 1, n
            deri(:, :, r, s) = deri(:, :, r, s) + matmul(eri_mo(:, :, r, s), u)
         end do
      end do
      do s = 1, n
         deri(:, :, :, s) = deri(:, :, :, s) &
                            + reshape(matmul(reshape(eri_mo(:, :, :, s), [n*n, n]), u), [n, n, n])
      end do
      ! Fourth index: ([pqr], t) -> ([pqr], s).
      deri = deri + reshape(matmul(reshape(eri_mo, [n*n*n, n]), u), [n, n, n, n])
      deri = deri + erix1
   end subroutine mp2_perturbed_eri

   subroutine mp2_perturbed_t2(deri, df, t2, orbital_energies, n_frozen, n_occ, dt2)
      !! The first-order response of the MP2 amplitudes, `d_X t2`
      !!
      !! Closed form, pycc's `_perturbed_t2` -- the amplitudes are not
      !! iterative, so their response is a divide:
      !!
      !!     d_X t2_ijab = [ d_X<ij|ab> + d_X f_ac t2_ijcb + d_X f_bc t2_ijac
      !!                     - d_X f_ik t2_kjab - d_X f_jk t2_ikab ] / D_ijab
      !!
      !! with `D_ijab = eps_i + eps_j - eps_a - eps_b`. The Fock blocks are the
      !! full oo/vv matrices of `d_X f`: the diagonal alone recovers the
      !! `-t2 d_X D` of a canonical gauge, and the off-diagonal rows are the
      !! non-canonical coupling this gauge carries instead.
      real(dp), intent(in) :: deri(:, :, :, :)   !! From `mp2_perturbed_eri`
      real(dp), intent(in) :: df(:, :)           !! From `mp2_perturbed_fock`
      real(dp), intent(in) :: t2(:, :, :, :)     !! (i, j, a, b), active
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_frozen, n_occ
      real(dp), allocatable, intent(out) :: dt2(:, :, :, :)

      real(dp) :: acc
      integer :: n_o, n_v, i, j, a, b, c, k

      n_o = size(t2, 1)
      n_v = size(t2, 3)
      allocate (dt2(n_o, n_o, n_v, n_v))
      do b = 1, n_v
         do a = 1, n_v
            do j = 1, n_o
               do i = 1, n_o
                  acc = deri(n_frozen + i, n_frozen + j, n_occ + a, n_occ + b)
                  do c = 1, n_v
                     acc = acc + df(n_occ + a, n_occ + c)*t2(i, j, c, b) &
                           + df(n_occ + b, n_occ + c)*t2(i, j, a, c)
                  end do
                  do k = 1, n_o
                     acc = acc - df(n_frozen + i, n_frozen + k)*t2(k, j, a, b) &
                           - df(n_frozen + j, n_frozen + k)*t2(i, k, a, b)
                  end do
                  dt2(i, j, a, b) = acc &
                                    /(orbital_energies(n_frozen + i) + orbital_energies(n_frozen + j) &
                                      - orbital_energies(n_occ + a) - orbital_energies(n_occ + b))
               end do
            end do
         end do
      end do
   end subroutine mp2_perturbed_t2

   subroutine mp2_perturbed_onepdm(t2, dt2, n_frozen, n_occ, n_mo, dd)
      !! The response of the unrelaxed one-particle density, `d_X gamma`
      !!
      !! The unrelaxed-density expressions differentiated by the product rule
      !! (pycc `_perturbed_unrelaxed_densities`, spin-adapted):
      !!
      !!     d_X gamma_ij = -( ta_imef l2_jmef + t2_imef la_jmef )
      !!     d_X gamma_ab =  ( ta_mnbe l2_mnae + t2_mnbe la_mnae )
      !!
      !! with `l2 = 2 (2 t2 - t2^{ab<->ba})` and `la` the same functional of
      !! `ta = d_X t2`. The cumulant's response needs no routine of its own:
      !! it is `mp2_cumulant_2pdm` evaluated at `ta`.
      real(dp), intent(in) :: t2(:, :, :, :), dt2(:, :, :, :)
      integer, intent(in) :: n_frozen, n_occ, n_mo
      real(dp), allocatable, intent(out) :: dd(:, :)

      real(dp), allocatable :: l2(:, :, :, :), la(:, :, :, :)
      real(dp) :: acc
      integer :: n_o, n_v, i, j, a, b, m, n, e, f

      n_o = size(t2, 1)
      n_v = size(t2, 3)
      allocate (l2(n_o, n_o, n_v, n_v), la(n_o, n_o, n_v, n_v))
      do b = 1, n_v
         do a = 1, n_v
            do j = 1, n_o
               do i = 1, n_o
                  l2(i, j, a, b) = 2.0_dp*(2.0_dp*t2(i, j, a, b) - t2(i, j, b, a))
                  la(i, j, a, b) = 2.0_dp*(2.0_dp*dt2(i, j, a, b) - dt2(i, j, b, a))
               end do
            end do
         end do
      end do

      allocate (dd(n_mo, n_mo))
      dd = 0.0_dp
      do j = 1, n_o
         do i = 1, n_o
            acc = 0.0_dp
            do f = 1, n_v
               do e = 1, n_v
                  do m = 1, n_o
                     acc = acc + dt2(i, m, e, f)*l2(j, m, e, f) &
                           + t2(i, m, e, f)*la(j, m, e, f)
                  end do
               end do
            end do
            dd(n_frozen + i, n_frozen + j) = -acc
         end do
      end do
      do b = 1, n_v
         do a = 1, n_v
            acc = 0.0_dp
            do e = 1, n_v
               do n = 1, n_o
                  do m = 1, n_o
                     acc = acc + dt2(m, n, b, e)*l2(m, n, a, e) &
                           + t2(m, n, b, e)*la(m, n, a, e)
                  end do
               end do
            end do
            dd(n_occ + a, n_occ + b) = acc
         end do
      end do
      deallocate (l2, la)
   end subroutine mp2_perturbed_onepdm

   subroutine mp2_perturbed_lagrangian(df, deri, dl, eri_mo, l_mo, orbital_energies, &
                                       d, dd, gam, dgam, n_occ, dip)
      !! The first-order response of the generalized-Fock Lagrangian, `d_X I'`
      !!
      !! pycc's `_perturbed_lagrangian`:
      !!
      !!     d_X I'_pq = -1/2 [ (d_X f (D + D^T))_pq + eps_p (dD_pq + dD_qp)
      !!                        + delta_{q occ} ( dD_rs L_rpsq + D_rs dL_rpsq
      !!                                        + dD_rs L_rqsp + D_rs dL_rqsp )
      !!                        + 4 ( d_X<pr|st> Gam_qrst + <pr|st> dGam_qrst ) ]
      !!
      !! **The one-electron derivative is the full matrix product
      !! `d_X f (D + D^T)`, never the diagonal `d_X eps` stencil** -- the
      !! off-diagonal `d_X f` couples the ov block a relaxed density carries,
      !! and the fixed-basis checks that would catch the shortcut share its
      !! assumption and are blind to it. Evaluated at the unrelaxed density and
      !! its response this drives the perturbed Z-vector; at the relaxed pair it
      !! is `d_X I`.
      real(dp), intent(in) :: df(:, :)             !! `d_X f`, full
      real(dp), intent(in), target, contiguous :: deri(:, :, :, :)
         !! `d_X <pq|rs>`, full
      real(dp), intent(in), target, contiguous :: dl(:, :, :, :)
         !! `2 deri - deri^{s<->r}`
      real(dp), intent(in), target, contiguous :: eri_mo(:, :, :, :), l_mo(:, :, :, :)
      real(dp), intent(in) :: orbital_energies(:)
      real(dp), intent(in) :: d(:, :), dd(:, :)    !! Density and its response
      real(dp), intent(in), target, contiguous :: gam(:, :, :, :), dgam(:, :, :, :)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: dip(:, :)

      real(dp), pointer :: deri2(:, :), gam2(:, :), eri2(:, :), dgam2(:, :)
      real(dp), allocatable :: two_e(:, :), dsym(:, :)
      real(dp) :: acc
      integer :: n_mo, n_cube, p, q, r, s

      n_mo = size(df, 1)
      n_cube = n_mo*n_mo*n_mo
      allocate (dip(n_mo, n_mo), two_e(n_mo, n_mo), dsym(n_mo, n_mo))

      ! The `4 (d_X<pr|st> Gam_qrst + <pr|st> dGam_qrst)` term is two matrix
      ! products, not a loop nest: `p` and `q` lead their arrays and Fortran
      ! stores columns first, so `(r,s,t)` is one contiguous trailing index and
      ! each contraction is `A B^T` over it. The remapping is a pointer view, so
      ! no `n_mo^4` copy is made to get it.
      deri2(1:n_mo, 1:n_cube) => deri
      gam2(1:n_mo, 1:n_cube) => gam
      eri2(1:n_mo, 1:n_cube) => eri_mo
      dgam2(1:n_mo, 1:n_cube) => dgam
      call pic_gemm(deri2, gam2, two_e, transb="T")
      call pic_gemm(eri2, dgam2, two_e, transb="T", beta=1.0_dp)

      ! `d_X f (D + D^T)` -- the full matrix product the docstring insists on.
      dsym = d + transpose(d)
      call pic_gemm(df, dsym, dip)

      do q = 1, n_mo
         do p = 1, n_mo
            dip(p, q) = dip(p, q) + orbital_energies(p)*(dd(p, q) + dd(q, p)) &
                        + 4.0_dp*two_e(p, q)
         end do
      end do

      ! This one stays a loop: `p` and `q` sit in the second and fourth slots of
      ! `l_mo`, so the `(r,s)` pair it contracts is strided and no reshape makes
      ! it a gemm without a transpose costing more than the contraction. It is
      ! `n_mo^4` against the term above's `n_mo^5`, so threading the independent
      ! `(p,q)` is the whole win.
      !$omp parallel do collapse(2) default(none) schedule(static) &
      !$omp    shared(dip, dd, d, l_mo, dl, n_mo, n_occ) private(p, q, r, s, acc)
      do q = 1, n_occ
         do p = 1, n_mo
            acc = 0.0_dp
            do s = 1, n_mo
               do r = 1, n_mo
                  acc = acc + dd(r, s)*(l_mo(r, p, s, q) + l_mo(r, q, s, p)) &
                        + d(r, s)*(dl(r, p, s, q) + dl(r, q, s, p))
               end do
            end do
            dip(p, q) = dip(p, q) + acc
         end do
      end do
      !$omp end parallel do
      dip = -0.5_dp*dip
   end subroutine mp2_perturbed_lagrangian

   subroutine mp2_perturbed_zvector_rhs(dip, df, dl, z, n_occ, rhs)
      !! The perturbed Z-vector right-hand side, `dX - G^X z`, solver-shaped
      !!
      !!     dX_ia  = d_X I'_ia - d_X I'_ai
      !!     (G^X z)_ia = (dL_ajib + dL_abij) z_jb + d_X f_ab z_ib - d_X f_ij z_ja
      !!
      !! `G^X z` is the perturbed orbital Hessian acting on the *unperturbed*
      !! Z-vector, which is what lets the perturbed solve reuse the unperturbed
      !! operator. Returned as `(n_vir, n_occ)` with the sign `cphf_solve`'s
      !! internal negation expects, so the response it returns is exactly the
      !! block **added** to `d_X D_rel`'s vo block.
      real(dp), intent(in) :: dip(:, :)      !! From `mp2_perturbed_lagrangian`
      real(dp), intent(in) :: df(:, :)
      real(dp), intent(in) :: dl(:, :, :, :)
      real(dp), intent(in) :: z(:, :)        !! Unperturbed Z-vector, (n_occ, n_vir)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: rhs(:, :)

      real(dp) :: acc
      integer :: n_mo, n_vir, i, a, j, b

      n_mo = size(dip, 1)
      n_vir = n_mo - n_occ
      allocate (rhs(n_vir, n_occ))
      do i = 1, n_occ
         do a = 1, n_vir
            acc = dip(i, n_occ + a) - dip(n_occ + a, i)
            do b = 1, n_vir
               do j = 1, n_occ
                  acc = acc - (dl(n_occ + a, j, i, n_occ + b) &
                               + dl(n_occ + a, n_occ + b, i, j))*z(j, b)
               end do
            end do
            do b = 1, n_vir
               acc = acc - df(n_occ + a, n_occ + b)*z(i, b)
            end do
            do j = 1, n_occ
               acc = acc + df(i, j)*z(j, a)
            end do
            rhs(a, i) = acc
         end do
      end do
   end subroutine mp2_perturbed_zvector_rhs

   subroutine mp2_perturbed_response(mol, coeff, orbital_energies, n_occ, &
                                     n_frozen, fx, sx, erix, mo1, t2, eri_mo, &
                                     l_mo, gam, drel, dt2, ddrel, di, error, &
                                     tol, zvec, xc, scf_density)
      !! The perturbed amplitudes, relaxed density and energy-weighted density
      !! for every nuclear perturbation
      !!
      !! pycc's `_perturbed_relaxed_density`, non-canonical, frozen-core
      !! aware. Per perturbation:
      !!
      !!     d_Y Drel = d_Y gamma + d_Y P_co - z^Y  (vo/ov from the solve)
      !!     G z^Y    = dX - G^Y z
      !!     d_Y I    = d_Y I'(Drel, d_Y Drel, Gam, d_Y Gam)
      !!
      !! **The core<->active derivative is a Sylvester relation, not a
      !! quotient rule** (Baeck, Watts and Bartlett, JCP 107, 3853 (1997),
      !! Eqs. 22-28). The gradient's divide
      !! `P_co[c,i] = (I'[c,i] - I'[i,c])/(eps_c - eps_i)` is the canonical
      !! collapse of `f P - P f = I' - I'^T`; a displacement leaves both
      !! blocks non-canonical, so the derivative keeps the off-diagonal
      !! `d_Y f` couplings:
      !!
      !!     d_Y P_co[c,i] = [ d_Y I'[c,i] - d_Y I'[i,c]
      !!                       - Sum_d d_Y f[c,d] P_co[d,i]
      !!                       + Sum_j P_co[c,j] d_Y f[j,i] ] / (eps_c - eps_i)
      !!
      !! `d` over the core, `j` over the active occupied. Dropping the
      !! off-diagonal terms is a real error, not a rounding one. The perturbed
      !! ov right-hand side picks up the derivative of the core<->active
      !! coupling as well: both the rotation's response against `L` and the
      !! unperturbed rotation against `d_Y L`.
      !!
      !! Two passes around one batched CPHF call. Pass 1 holds one
      !! perturbation's `nmo^4` derivative at a time and folds out what pass 2
      !! needs from it -- the Z-vector right-hand side and `d_Y I` at the
      !! *pre-solve* density response. Pass 2 adds the solved blocks; their
      !! share of `d_Y I` involves only the unperturbed `L`, because the
      !! response enters the Lagrangian linearly, which is what lets the tensor
      !! be dropped before the solve rather than rebuilt after it.
      !!
      !! The cumulant's response is not stacked: it is `mp2_cumulant_2pdm`
      !! evaluated at `dt2(:, :, :, :, y)`, and `o^2 v^2` amplitudes are the
      !! compact carrier where `nmo^4` cumulants are not.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ
      integer, intent(in) :: n_frozen
      real(dp), intent(in) :: fx(:, :, :), sx(:, :, :)   !! From `mp2_first_order_skeletons`
      real(dp), intent(in) :: erix(:, :, :, :, :)
      real(dp), intent(in) :: mo1(:, :, :, :)   !! (n_mo, n_occ, 3, natm)
      real(dp), intent(in) :: t2(:, :, :, :)    !! Active-occupied amplitudes
      real(dp), intent(in) :: eri_mo(:, :, :, :), l_mo(:, :, :, :)
      real(dp), intent(in) :: gam(:, :, :, :)   !! From `mp2_cumulant_2pdm`
      real(dp), intent(in) :: drel(:, :)        !! Relaxed 1-PDM, MO
      real(dp), allocatable, intent(out) :: dt2(:, :, :, :, :)
      real(dp), allocatable, intent(out) :: ddrel(:, :, :)
      real(dp), allocatable, intent(out) :: di(:, :, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tol
         !! Solver tolerance for both Z-vector solves. Defaults to 1e-13,
         !! because the references it is gated against are dense solves and the
         !! whole iterative floor in a comparison is this side's.
      real(dp), allocatable, intent(out), optional :: zvec(:, :)
         !! The unperturbed Z-vector `z_ia`, `(n_occ, n_vir)`, for a caller
         !! that wants to cross-check it against the gradient's relaxed
         !! density.
      type(xc_context_t), intent(inout), optional :: xc
         !! The functional the reference used. Present, the two Z-vector
         !! solves become Kohn-Sham ones, the derivative weight `dl` carries
         !! the exchange fraction, and every `L`-fold's kernel share is added
         !! through the operator applies at the bottom of this file -- the
         !! full-derivative folds through `ks_dl_pair_kernel`, whose `g_xc`
         !! term is where the kernel's own density dependence enters. `l_mo`
         !! must arrive already `k`-scaled.
      real(dp), intent(in), optional :: scf_density(:, :)
         !! The converged reference density. Required with `xc`.

      real(dp), allocatable :: d0(:, :), ip0(:, :), rhs0(:, :, :), zsol(:, :, :)
      real(dp), allocatable :: z(:, :), u(:, :), df(:, :), deri(:, :, :, :)
      real(dp), allocatable :: dl(:, :, :, :), ta(:, :, :, :), ddg(:, :)
      real(dp), allocatable :: dgam_y(:, :, :, :), dip(:, :), rhs1(:, :)
      real(dp), allocatable :: rhs(:, :, :), zx(:, :, :), zm(:, :)
      real(dp), allocatable :: pco(:, :), dpco(:, :)
      real(dp), allocatable :: vk_d0(:, :), vk_drel(:, :), vk_z(:, :), vk_pco(:, :)
      real(dp), allocatable :: kd_d0(:, :, :), kd_drel(:, :, :), kd_z(:, :, :)
      real(dp), allocatable :: kd_pco(:, :, :)
      real(dp), allocatable :: z_emb(:, :), pco_emb(:, :), mx(:, :), vk_mx(:, :)
      real(dp), allocatable :: ddmat_ao(:, :), w_d0(:, :), w_drel(:, :), w_z(:, :)
      real(dp), allocatable :: w_pco(:, :), vk_dd(:, :), vk_dpco(:, :), vk_zm(:, :)
      real(dp), allocatable :: xc_occ(:, :), dpco_emb(:, :)
      logical :: dh
      real(dp) :: use_tol, acc, kf
      integer :: n_ao, n_mo, n_vir, n_act, n_pert, y, atom, comp, i, j, a, b, c, d
      integer :: m, e, f, p, q, r, s

      if (error%has_error()) return

      dh = present(xc)
      kf = 1.0_dp
      if (dh) then
         if (.not. present(scf_density)) then
            call error%set(ERROR_VALIDATION, "the double-hybrid perturbed response "// &
                           "needs the converged reference density: the kernel and "// &
                           "its derivatives are evaluated at it")
            return
         end if
         kf = xc%exx_fraction
      end if

      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)
      n_vir = n_mo - n_occ
      n_act = n_occ - n_frozen
      n_pert = size(sx, 3)
      if (size(t2, 1) /= n_act) then
         call error%set(ERROR_VALIDATION, "the perturbed response was handed "// &
                        "amplitudes whose occupied count does not match "// &
                        "n_occ - n_frozen.")
         return
      end if
      use_tol = 1.0e-13_dp
      if (present(tol)) use_tol = tol

      ! ---- the unperturbed pieces the response reuses ----------------------
      ! The unrelaxed density from the amplitudes (active blocks only -- the
      ! core rows carry no correlation), its Lagrangian, the core<->active
      ! divide, and the Z-vector: z = G^-1 X, X_ia = I'_ia - I'_ai over the full
      ! occupied space, minus the core-rotation coupling when a core exists.
      allocate (d0(n_mo, n_mo))
      d0 = 0.0_dp
      do j = 1, n_act
         do i = 1, n_act
            acc = 0.0_dp
            do f = 1, n_vir
               do e = 1, n_vir
                  do m = 1, n_act
                     acc = acc + t2(i, m, e, f)*(2.0_dp*(2.0_dp*t2(j, m, e, f) - t2(j, m, f, e)))
                  end do
               end do
            end do
            d0(n_frozen + i, n_frozen + j) = -acc
         end do
      end do
      do b = 1, n_vir
         do a = 1, n_vir
            acc = 0.0_dp
            do e = 1, n_vir
               do j = 1, n_act
                  do m = 1, n_act
                     acc = acc + t2(m, j, b, e)*(2.0_dp*(2.0_dp*t2(m, j, a, e) - t2(m, j, e, a)))
                  end do
               end do
            end do
            d0(n_occ + a, n_occ + b) = acc
         end do
      end do

      call mp2_mo_lagrangian(eri_mo, orbital_energies, d0, gam, n_occ, ip0)
      if (dh) then
         ! The unrelaxed Lagrangian's Kohn-Sham lift, outside the routine for
         ! the module header's bit-for-bit reason.
         allocate (vk_d0(n_mo, n_mo))
         call ks_pair_kernel(xc, mol, coeff, scf_density, d0, vk_d0, error)
         if (error%has_error()) return
         call ks_occ_fold_patch(eri_mo, d0, n_occ, kf, vk_d0, ip0)
      end if
      if (n_frozen > 0) then
         allocate (pco(n_frozen, n_act))
         do i = 1, n_act
            do c = 1, n_frozen
               pco(c, i) = (ip0(c, n_frozen + i) - ip0(n_frozen + i, c)) &
                           /(orbital_energies(c) - orbital_energies(n_frozen + i))
            end do
         end do
         if (dh) then
            ! The rotation as the pair folds see it: one triangle of the
            ! core<->active divide, the operator applies symmetrising it.
            allocate (pco_emb(n_mo, n_mo), vk_pco(n_mo, n_mo))
            pco_emb = 0.0_dp
            do i = 1, n_act
               do c = 1, n_frozen
                  pco_emb(n_frozen + i, c) = pco(c, i)
               end do
            end do
            call ks_pair_kernel(xc, mol, coeff, scf_density, pco_emb, vk_pco, error)
            if (error%has_error()) return
         end if
      end if
      allocate (rhs0(n_vir, n_occ, 1))
      do i = 1, n_occ
         do a = 1, n_vir
            acc = ip0(i, n_occ + a) - ip0(n_occ + a, i)
            if (n_frozen > 0) then
               ! The core rotation rides the ov right-hand side through the
               ! orbital Hessian: z_jc = -P_co^T against the unperturbed L.
               do c = 1, n_frozen
                  do j = 1, n_act
                     acc = acc + pco(c, j)*(l_mo(n_occ + a, n_frozen + j, i, c) &
                                            + l_mo(n_occ + a, c, i, n_frozen + j))
                  end do
               end do
               if (dh) acc = acc + vk_pco(n_occ + a, i)
            end if
            rhs0(a, i, 1) = -acc
         end do
      end do
      ! `xc` and `scf_density` forward as they stand: absent, `cphf_solve`
      ! sees the Hartree-Fock equations it always solved.
      call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zsol, &
                      error=error, mo_rhs=rhs0, tol=use_tol, max_iter=200, &
                      xc=xc, density=scf_density)
      if (error%has_error()) return
      allocate (z(n_occ, n_vir))
      z = transpose(zsol(:, :, 1))
      deallocate (rhs0, zsol, ip0)
      if (present(zvec)) zvec = z

      ! ---- the kernel-derivative stacks, one grid pass per density ---------
      !
      ! Every `d_X L` fold below contracts one of four fixed densities, and
      ! the skeleton share of its kernel term depends on the density alone,
      ! so it is built here for all `3 natm` perturbations at once. The
      ! per-perturbation shares -- the rotations and the `g_xc` term -- are
      ! assembled inside the loop, where `U^X` first exists.
      if (dh) then
         allocate (z_emb(n_mo, n_mo), vk_z(n_mo, n_mo), vk_drel(n_mo, n_mo))
         z_emb = 0.0_dp
         z_emb(1:n_occ, n_occ + 1:n_mo) = z
         call ks_pair_kernel(xc, mol, coeff, scf_density, z_emb, vk_z, error)
         call ks_pair_kernel(xc, mol, coeff, scf_density, drel, vk_drel, error)
         call ks_kernel_deriv_stack(xc, mol, coeff, scf_density, d0, kd_d0, error)
         call ks_kernel_deriv_stack(xc, mol, coeff, scf_density, drel, kd_drel, error)
         call ks_kernel_deriv_stack(xc, mol, coeff, scf_density, z_emb, kd_z, error)
         if (n_frozen > 0) then
            call ks_kernel_deriv_stack(xc, mol, coeff, scf_density, pco_emb, &
                                       kd_pco, error)
         end if
         if (error%has_error()) return
         allocate (mx(n_mo, n_mo), vk_mx(n_mo, n_mo), ddmat_ao(n_ao, n_ao))
         allocate (w_d0(n_mo, n_mo), w_drel(n_mo, n_mo), w_z(n_mo, n_mo))
         allocate (vk_dd(n_mo, n_mo), vk_zm(n_mo, n_mo), xc_occ(n_mo, n_mo))
         if (n_frozen > 0) then
            allocate (w_pco(n_mo, n_mo), vk_dpco(n_mo, n_mo))
            allocate (dpco_emb(n_mo, n_mo))
         end if
      end if

      ! ---- pass 1: one perturbation's derivative resident at a time --------
      allocate (dt2(n_act, n_act, n_vir, n_vir, n_pert))
      allocate (ddrel(n_mo, n_mo, n_pert), di(n_mo, n_mo, n_pert))
      allocate (rhs(n_vir, n_occ, n_pert))
      do y = 1, n_pert
         atom = (y - 1)/3 + 1
         comp = y - 3*(atom - 1)
         if (dh) then
            ! The kernel's response to this perturbation's first-order
            ! density: `mo1` is the carrier -- `-S^(X)/2` occupied rows, the
            ! solved response below -- and the Brillouin rewrite inside
            ! `mp2_full_u` redistributes only the antisymmetric part of the
            ! core<->active block, which the symmetrised apply never sees.
            ! One matrix therefore serves the rewrite and `d_X f` alike.
            mx = 0.0_dp
            mx(:, 1:n_occ) = mo1(:, :, comp, atom)
            call ks_pair_kernel(xc, mol, coeff, scf_density, mx, vk_mx, error)
            if (error%has_error()) return
         end if
         call mp2_full_u(mo1(:, :, comp, atom), sx(:, :, y), n_occ, u, &
                         n_frozen=n_frozen, fx1=fx(:, :, y), l_mo=l_mo, &
                         orbital_energies=orbital_energies)
         if (dh .and. n_frozen > 0) then
            ! The Brillouin rewrite's kernel share, folded into the rotation
            ! after the fact: the rewrite's divides are independent -- each
            ! reads only the solved virtual rows -- so lifting `acc` by
            ! `vk_mx(c, i)` inside and correcting `u` here are the same
            ! rotation, and the routine's text stays the Hartree-Fock one
            ! (module header).
            do i = n_frozen + 1, n_occ
               do c = 1, n_frozen
                  u(c, i) = u(c, i) - vk_mx(c, i) &
                            /(orbital_energies(c) - orbital_energies(i))
                  u(i, c) = -sx(c, i, y) - u(c, i)
               end do
            end do
         end if
         call mp2_perturbed_fock(fx(:, :, y), sx(:, :, y), u, l_mo, &
                                 orbital_energies, n_occ, df)
         ! The mean-field response's kernel share: every element, since the
         ! perturbed Fock's folds run over the full square.
         if (dh) df = df + vk_mx
         call mp2_perturbed_eri(erix(:, :, :, :, y), u, eri_mo, deri)
         ! Two copies of the weight build, split on the exchange fraction,
         ! for the reason `mp2_first_order_skeletons` gives.
         allocate (dl(n_mo, n_mo, n_mo, n_mo))
         if (kf == 1.0_dp) then
            do s = 1, n_mo
               do r = 1, n_mo
                  dl(:, :, r, s) = 2.0_dp*deri(:, :, r, s) - deri(:, :, s, r)
               end do
            end do
         else
            do s = 1, n_mo
               do r = 1, n_mo
                  dl(:, :, r, s) = 2.0_dp*deri(:, :, r, s) - kf*deri(:, :, s, r)
               end do
            end do
         end if
         if (dh) then
            ! The full-derivative kernel shares for this perturbation. The
            ! reference density's matrix response is `4 C sym(mx) C^T` --
            ! `nuclear_apply`'s own first-order density -- and feeds the
            ! `g_xc` term of every `d_X L` fold.
            ddmat_ao = 2.0_dp*matmul(coeff, &
                                     matmul(mx + transpose(mx), transpose(coeff)))
            call ks_dl_pair_kernel(xc, mol, coeff, scf_density, d0, vk_d0, &
                                   kd_d0(:, :, y), u, ddmat_ao, w_d0, error)
            call ks_dl_pair_kernel(xc, mol, coeff, scf_density, drel, vk_drel, &
                                   kd_drel(:, :, y), u, ddmat_ao, w_drel, error)
            call ks_dl_pair_kernel(xc, mol, coeff, scf_density, z_emb, vk_z, &
                                   kd_z(:, :, y), u, ddmat_ao, w_z, error)
            if (n_frozen > 0) then
               call ks_dl_pair_kernel(xc, mol, coeff, scf_density, pco_emb, &
                                      vk_pco, kd_pco(:, :, y), u, ddmat_ao, &
                                      w_pco, error)
            end if
            if (error%has_error()) return
         end if
         call mp2_perturbed_t2(deri, df, t2, orbital_energies, n_frozen, n_occ, ta)
         dt2(:, :, :, :, y) = ta
         call mp2_perturbed_onepdm(t2, ta, n_frozen, n_occ, n_mo, ddg)
         call mp2_cumulant_2pdm(ta, n_frozen, n_occ, n_mo, dgam_y)
         ! The Z-vector right-hand side wants the Lagrangian's response at the
         ! unrelaxed pair; the energy-weighted density wants it at the relaxed
         ! one. Both while the derivative tensor is still here.
         call mp2_perturbed_lagrangian(df, deri, dl, eri_mo, l_mo, &
                                       orbital_energies, d0, ddg, gam, dgam_y, &
                                       n_occ, dip)
         if (dh) then
            ! The occupied-column folds' kernel shares -- `dD . L` through the
            ! apply at `ddg`, `D . d_X L` through the full-derivative helper
            ! -- added with the routine's own `-1/2` after the fact, so its
            ! Hartree-Fock text stands (module header).
            call ks_pair_kernel(xc, mol, coeff, scf_density, ddg, vk_dd, error)
            if (error%has_error()) return
            xc_occ = vk_dd + w_d0
            dip(:, 1:n_occ) = dip(:, 1:n_occ) - 0.5_dp*xc_occ(:, 1:n_occ)
         end if
         call mp2_perturbed_zvector_rhs(dip, df, dl, z, n_occ, rhs1)
         ! The `dL . z` fold's kernel share, subtracted the way the explicit
         ! fold is.
         if (dh) rhs1 = rhs1 - w_z(n_occ + 1:n_mo, 1:n_occ)
         if (n_frozen > 0) then
            ! The Sylvester derivative of the core<->active divide (header):
            ! the off-diagonal d_Y f couplings within the core and active
            ! blocks, then the rotation's response folded into the ov right-hand
            ! side alongside the unperturbed rotation against d_Y L.
            allocate (dpco(n_frozen, n_act))
            do i = 1, n_act
               do c = 1, n_frozen
                  acc = dip(c, n_frozen + i) - dip(n_frozen + i, c)
                  do d = 1, n_frozen
                     acc = acc - df(c, d)*pco(d, i)
                  end do
                  do j = 1, n_act
                     acc = acc + pco(c, j)*df(n_frozen + j, n_frozen + i)
                  end do
                  dpco(c, i) = acc/(orbital_energies(c) - orbital_energies(n_frozen + i))
               end do
            end do
            if (dh) then
               dpco_emb = 0.0_dp
               do i = 1, n_act
                  do c = 1, n_frozen
                     dpco_emb(n_frozen + i, c) = dpco(c, i)
                  end do
               end do
               call ks_pair_kernel(xc, mol, coeff, scf_density, dpco_emb, &
                                   vk_dpco, error)
               if (error%has_error()) return
            end if
            do i = 1, n_occ
               do a = 1, n_vir
                  acc = 0.0_dp
                  do c = 1, n_frozen
                     do j = 1, n_act
                        acc = acc + dpco(c, j)*(l_mo(n_occ + a, n_frozen + j, i, c) &
                                                + l_mo(n_occ + a, c, i, n_frozen + j)) &
                              + pco(c, j)*(dl(n_occ + a, n_frozen + j, i, c) &
                                           + dl(n_occ + a, c, i, n_frozen + j))
                     end do
                  end do
                  if (dh) acc = acc + vk_dpco(n_occ + a, i) + w_pco(n_occ + a, i)
                  rhs1(a, i) = rhs1(a, i) + acc
               end do
            end do
            do i = 1, n_act
               do c = 1, n_frozen
                  ddg(c, n_frozen + i) = ddg(c, n_frozen + i) + dpco(c, i)
                  ddg(n_frozen + i, c) = ddg(n_frozen + i, c) + dpco(c, i)
               end do
            end do
            deallocate (dpco)
         end if
         rhs(:, :, y) = rhs1
         deallocate (dip)
         call mp2_perturbed_lagrangian(df, deri, dl, eri_mo, l_mo, &
                                       orbital_energies, drel, ddg, gam, &
                                       dgam_y, n_occ, dip)
         if (dh) then
            ! As at the unrelaxed pair above; `ddg` may have just taken the
            ! core<->active response, so its kernel share is re-applied at
            ! the matrix the fold read.
            if (n_frozen > 0) then
               call ks_pair_kernel(xc, mol, coeff, scf_density, ddg, vk_dd, error)
               if (error%has_error()) return
            end if
            xc_occ = vk_dd + w_drel
            dip(:, 1:n_occ) = dip(:, 1:n_occ) - 0.5_dp*xc_occ(:, 1:n_occ)
         end if
         di(:, :, y) = dip
         ddrel(:, :, y) = ddg
         deallocate (u, df, deri, dl, ta, ddg, dgam_y, dip, rhs1)
      end do

      call cphf_solve(mol, coeff, orbital_energies, n_occ, response=zx, &
                      error=error, mo_rhs=rhs, tol=use_tol, max_iter=200, &
                      xc=xc, density=scf_density)
      if (error%has_error()) return
      deallocate (rhs)

      ! ---- pass 2: the solved blocks, and their linear share of d_Y I ------
      allocate (zm(n_mo, n_mo))
      do y = 1, n_pert
         ddrel(n_occ + 1:n_mo, 1:n_occ, y) = ddrel(n_occ + 1:n_mo, 1:n_occ, y) + zx(:, :, y)
         ddrel(1:n_occ, n_occ + 1:n_mo, y) = ddrel(1:n_occ, n_occ + 1:n_mo, y) &
                                             + transpose(zx(:, :, y))
         zm = 0.0_dp
         zm(n_occ + 1:n_mo, 1:n_occ) = zx(:, :, y)
         zm(1:n_occ, n_occ + 1:n_mo) = transpose(zx(:, :, y))
         do q = 1, n_mo
            do p = 1, n_mo
               di(p, q, y) = di(p, q, y) &
                             - 0.5_dp*orbital_energies(p)*(zm(p, q) + zm(q, p))
            end do
         end do
         do q = 1, n_occ
            do p = 1, n_mo
               acc = 0.0_dp
               do s = 1, n_mo
                  do r = 1, n_mo
                     acc = acc + zm(r, s)*(l_mo(r, p, s, q) + l_mo(r, q, s, p))
                  end do
               end do
               di(p, q, y) = di(p, q, y) - 0.5_dp*acc
            end do
         end do
         if (dh) then
            ! The solved blocks' own kernel share of `d_Y I`, unperturbed
            ! operator only -- the same linearity that lets pass 2 run after
            ! the derivative tensors are gone.
            call ks_pair_kernel(xc, mol, coeff, scf_density, zm, vk_zm, error)
            if (error%has_error()) return
            di(:, 1:n_occ, y) = di(:, 1:n_occ, y) - 0.5_dp*vk_zm(:, 1:n_occ)
         end if
      end do
      deallocate (zm, zx, d0, z)
      if (allocated(pco)) deallocate (pco)
      if (allocated(vk_d0)) then
         deallocate (vk_d0, vk_drel, vk_z, kd_d0, kd_drel, kd_z)
         deallocate (z_emb, mx, vk_mx, ddmat_ao, w_d0, w_drel, w_z)
         deallocate (vk_dd, vk_zm, xc_occ)
      end if
      if (allocated(vk_pco)) deallocate (vk_pco, pco_emb, kd_pco, w_pco, &
                                         vk_dpco, dpco_emb)
   end subroutine mp2_perturbed_response

   subroutine mp2_correlation_hessian(mol, coeff, orbital_energies, density, &
                                      n_occ, n_frozen, hess_corr, hess_ref, &
                                      error, tol, xc, pt2_scale, dipole_derivatives)
      !! The analytic MP2 correlation Hessian, `(3, 3, natm, natm)`
      !!
      !! The three groups of pycc's `_hessian_blocks` (`route == 'aod'`)
      !! combined: the fixed-density second skeletons (pass 1,
      !! `mp2_skeleton_hessian`), the orbital response
      !! (`mp2_orbital_response_term`), and the 2n+1 density response
      !!
      !!     d_Y D_rel . f^(X) + d_Y Gamma . <pq|rs>^(X) + d_Y I . S^(X)
      !!
      !! from `mp2_perturbed_response` (pass 2). **Two passes, not one fused
      !! loop**: the second-derivative block never meets `erix`/`d_Y Gamma` in a
      !! contraction, so the `nao^4` effective density is built, swept and freed
      !! before any pass-2 tensor is allocated -- peak `max(9, 6+6)` working
      !! sets instead of `9+6+6`.
      !!
      !! `hess_corr` is the correlation block alone. `hess_ref` is the SCF
      !! reference's AO-dependent skeleton, deposited from pass 1's shared
      !! integral sweep; the caller completes the reference by delegating its
      !! CPHF response to `response_hessian`, so the total Hessian is
      !!
      !!     nuclear_repulsion_hessian + hess_ref + response_hessian + hess_corr
      !!
      !! Frozen-core aware: `n_frozen > 0` runs the core<->active rewrite of
      !! `U^X`, the pair augmentation and the Sylvester derivative;
      !! `n_frozen = 0` reproduces the all-electron path verbatim.
      !!
      !! **With `xc`, this returns something else**: the correlation block of a
      !! double hybrid, over the Kohn-Sham reference whose orbitals and
      !! `density` arrive here, scaled by `pt2_scale`. `hess_ref` is then still
      !! the Hartree-Fock-shaped skeleton and is to be **discarded**: a
      !! double-hybrid caller takes its whole reference block from `ks_hessian`
      !! instead.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)            !! C, (n_ao, n_mo)
      real(dp), intent(in) :: orbital_energies(:)    !! (n_mo), Hartree
      real(dp), intent(in) :: density(:, :)          !! Converged AO density
      integer, intent(in) :: n_occ                   !! Doubly occupied count
      integer, intent(in) :: n_frozen                !! Frozen-core count
      real(dp), allocatable, intent(out) :: hess_corr(:, :, :, :)
      real(dp), allocatable, intent(out) :: hess_ref(:, :, :, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: tol
         !! Solver tolerance for the CPHF and both Z-vector solves. Defaults to
         !! 1e-13, below the solvers' own defaults, because the references this
         !! assembly is gated against are dense solves.
      type(xc_context_t), intent(inout), optional :: xc
         !! The functional the Kohn-Sham reference used. Present switches the
         !! whole ladder to the double hybrid's correlation block -- see the
         !! header. `density` is then where every kernel is evaluated.
      real(dp), intent(in), optional :: pt2_scale
         !! The functional's PT2 coefficient. Applied once, to the finished
         !! correlation block: the perturbed Z-vector solves are linear, so
         !! scaling their inputs and scaling their outputs are the same
         !! operation and doing both is the error.
      real(dp), allocatable, intent(out), optional :: dipole_derivatives(:, :)
         !! The perturbative term's share of `d mu_a / dR_(X,b)`, `(3, 3*natm)`,
         !! scaled by `pt2_scale` like the correlation block itself. A double
         !! hybrid's dipole is the relaxed one, so its infrared intensity needs
         !! this on top of the reference's derivative, exactly as its Hessian
         !! needs `hess_corr` on top of `ks_hessian`.
         !!
         !! Requires `xc`, and errors without it: over a Hartree-Fock reference
         !! an MP2 dipole derivative is a different quantity with a different
         !! assembly.

      real(dp), allocatable :: gradient(:, :), dm1mo(:, :), w_ao(:, :)
      real(dp), allocatable :: eri_packed(:, :), ovov(:, :, :, :), t2(:, :, :, :)
      real(dp), allocatable :: gamma_eff(:, :, :, :)
      real(dp), allocatable :: fx(:, :, :), sx(:, :, :), erix(:, :, :, :, :)
      real(dp), allocatable :: eri_mo(:, :, :, :), l_mo(:, :, :, :)
      real(dp), allocatable :: gam(:, :, :, :), imat(:, :)
      real(dp), allocatable :: ip(:, :), xov(:, :), i2(:, :)
      real(dp), allocatable :: xt(:, :), it(:, :), pf(:, :)
      real(dp), allocatable :: xov_st(:, :, :), i2_st(:, :, :), pf_st(:, :, :)
      real(dp), allocatable :: eip1(:, :, :, :, :), h1(:, :, :, :), s1(:, :, :, :)
      real(dp), allocatable :: h1a(:, :, :), s1a(:, :, :), mo1(:, :, :, :)
      real(dp), allocatable :: orb(:, :), dt2(:, :, :, :, :)
      real(dp), allocatable :: ddrel(:, :, :), di(:, :, :), dgam(:, :, :, :)
      real(dp), allocatable :: vk_drel(:, :), kd_drel_mo(:, :, :), vkpf(:, :)
      logical :: dh
      real(dp) :: use_tol, resp, kf, cscale
      integer :: n_ao, n_mo, natm, n_pert, a, x, ix, iy, ax, by, cx, cy
      integer :: p, q, r, s

      if (error%has_error()) return
      use_tol = 1.0e-13_dp
      if (present(tol)) use_tol = tol

      ! The double hybrid's settings, and the ordinary MP2 values that leave
      ! every expression below bit for bit what it was: full exchange, unit
      ! scale, no kernel anywhere.
      dh = present(xc)
      kf = 1.0_dp
      cscale = 1.0_dp
      if (dh) then
         kf = xc%exx_fraction
         if (present(pt2_scale)) cscale = pt2_scale
      end if

      n_ao = mol%nao
      n_mo = size(coeff, 2)
      natm = mol%natm
      n_pert = 3*natm

      ! ---- the unperturbed ladder, from the gradient ------------------------
      if (dh) then
         call czt_mp2_gradient(mol, coeff, orbital_energies, n_occ, gradient, &
                               error, n_frozen=n_frozen, &
                               relaxed_density_mo=dm1mo, &
                               energy_weighted_ao=w_ao, xc=xc, &
                               scf_density=density, pt2_scale=cscale)
      else
         call czt_mp2_gradient(mol, coeff, orbital_energies, n_occ, gradient, &
                               error, n_frozen=n_frozen, &
                               relaxed_density_mo=dm1mo, &
                               energy_weighted_ao=w_ao)
      end if
      if (error%has_error()) return

      call mol%eris_packed(eri_packed)
      call transform_ovov(eri_packed, coeff, n_frozen, n_occ, n_ao, n_mo, ovov)
      call build_amplitudes(ovov, orbital_energies, n_frozen, n_occ - n_frozen, &
                            n_mo - n_occ, n_occ, t2)
      deallocate (eri_packed, ovov)

      ! ---- pass 1: the fixed-density second skeletons, both blocks ----------
      call build_effective_2pdm_ao(t2, dm1mo, coeff, n_ao, n_mo, n_occ, n_frozen, &
                                   gamma_eff, k_scale=kf)
      ! `xc` forwards as it stands -- absent, the sweep is the Hartree-Fock
      ! one and `scf_density` goes unread.
      call mp2_skeleton_hessian(mol, gamma_eff, dm1mo, w_ao, coeff, &
                                orbital_energies, n_occ, hess_corr, hess_ref, &
                                error, xc=xc, scf_density=density)
      if (error%has_error()) return
      ! Freed before any pass-2 tensor exists -- the two-pass memory shape.
      deallocate (gamma_eff)

      ! ---- pass 2 inputs: skeletons, carriers, and the CPHF ----------------
      call mp2_first_order_skeletons(mol, coeff, n_occ, fx, sx, erix, error)
      if (error%has_error()) return
      ! The Kohn-Sham lift lives outside the skeleton builder on purpose:
      ! the Hartree-Fock routine's text is part of the bit-for-bit contract
      ! (module header), and everything the lift needs -- the occupied trace
      ! of `erix` and the potential's grid derivative -- is in hand here.
      if (dh) then
         call mp2_ks_first_order_fock(mol, xc, density, coeff, n_occ, erix, kf, &
                                      fx, error)
         if (error%has_error()) return
      end if
      call mp2_mo_eri_physicist(mol, coeff, eri_mo, error)
      if (error%has_error()) return

      ! Two copies of the weight build, split on the exchange fraction, for
      ! the reason `mp2_first_order_skeletons` gives: full exchange keeps its
      ! literal statements bit for bit.
      allocate (l_mo(n_mo, n_mo, n_mo, n_mo))
      if (kf == 1.0_dp) then
         do s = 1, n_mo
            do r = 1, n_mo
               do q = 1, n_mo
                  do p = 1, n_mo
                     l_mo(p, q, r, s) = 2.0_dp*eri_mo(p, q, r, s) - eri_mo(p, q, s, r)
                  end do
               end do
            end do
         end do
      else
         do s = 1, n_mo
            do r = 1, n_mo
               do q = 1, n_mo
                  do p = 1, n_mo
                     l_mo(p, q, r, s) = 2.0_dp*eri_mo(p, q, r, s) &
                                        - kf*eri_mo(p, q, s, r)
                  end do
               end do
            end do
         end do
      end if

      call mp2_cumulant_2pdm(t2, n_frozen, n_occ, n_mo, gam)
      call mp2_mo_lagrangian(eri_mo, orbital_energies, dm1mo, gam, n_occ, imat)
      if (dh) then
         ! The Lagrangian's occupied-column fold, lifted to the reference's
         ! exchange fraction plus its kernel -- outside `mp2_mo_lagrangian`
         ! for the bit-for-bit reason above.
         allocate (vk_drel(n_mo, n_mo))
         call ks_pair_kernel(xc, mol, coeff, density, dm1mo, vk_drel, error)
         if (error%has_error()) return
         call ks_occ_fold_patch(eri_mo, dm1mo, n_occ, kf, vk_drel, imat)
      end if

      ! All-electron non-canonical, the augmentation is the exact identity, so
      ! `mp2_pair_rotation_augment` is skipped. A frozen core forces it: the
      ! independent core<->active rotation has no CPHF counterpart and its
      ! `P^(X)` carrier must meet `f^(Y)` in the orbital response.
      allocate (xov_st(n_mo, n_mo, n_pert), i2_st(n_mo, n_mo, n_pert))
      if (n_frozen > 0) allocate (pf_st(n_mo, n_mo, n_pert))
      if (dh) then
         ! The skeleton derivative of the Lagrangian's kernel fold, all
         ! perturbations from one grid pass.
         call ks_kernel_deriv_stack(xc, mol, coeff, density, dm1mo, &
                                    kd_drel_mo, error)
         if (error%has_error()) return
         if (n_frozen > 0) allocate (vkpf(n_mo, n_mo))
      end if
      do x = 1, n_pert
         call mp2_skeleton_lagrangian(fx(:, :, x), sx(:, :, x), &
                                      erix(:, :, :, :, x), dm1mo, gam, imat, &
                                      n_occ, ip, xov, i2)
         if (dh) then
            ! The fold's Kohn-Sham lift -- exchange fraction plus the
            ! kernel's skeleton derivative -- lands on `ip`, and the two
            ! carriers are rebuilt from it exactly as the routine builds
            ! them. Note `fx` and `imat` already arrived lifted, so those
            ! terms of `I'^(X)` need nothing here.
            call ks_occ_fold_patch(erix(:, :, :, :, x), dm1mo, n_occ, kf, &
                                   kd_drel_mo(:, :, x), ip)
            xov = transpose(ip) - ip
            i2 = ip
            i2(n_occ + 1:n_mo, 1:n_occ) = transpose(ip(1:n_occ, n_occ + 1:n_mo))
         end if
         if (n_frozen > 0) then
            call mp2_pair_rotation_augment(ip, xov, i2, l_mo, orbital_energies, &
                                           n_occ, n_frozen, .false., xt, it, pf)
            if (dh) then
               ! The kernel share of the augmentation's pair sums -- the same
               ! `P^(X)` against `L` folds, as one operator apply, patched
               ! where the augmented carriers were just assembled so the
               ! transcription above stays checkable against pycc line by
               ! line.
               call ks_pair_kernel(xc, mol, coeff, density, pf, vkpf, error)
               if (error%has_error()) return
               xt(n_occ + 1:n_mo, 1:n_occ) = xt(n_occ + 1:n_mo, 1:n_occ) &
                                             + 0.5_dp*vkpf(n_occ + 1:n_mo, 1:n_occ)
               it(1:n_occ, 1:n_occ) = it(1:n_occ, 1:n_occ) &
                                      - 0.5_dp*vkpf(1:n_occ, 1:n_occ)
            end if
            xov_st(:, :, x) = xt
            i2_st(:, :, x) = it
            pf_st(:, :, x) = pf
            deallocate (xt, it, pf)
         else
            xov_st(:, :, x) = xov
            i2_st(:, :, x) = i2
         end if
         deallocate (ip, xov, i2)
      end do

      call eri_ip1_block(mol, eip1, error)
      if (error%has_error()) return
      allocate (h1(n_ao, n_ao, 3, natm), s1(n_ao, n_ao, 3, natm))
      do a = 1, natm
         call make_h1_atom(mol, density, eip1, a, h1a, error, k_scale=kf)
         call overlap_deriv_atom(mol, a, s1a, error)
         if (error%has_error()) return
         h1(:, :, :, a) = h1a
         s1(:, :, :, a) = s1a
         deallocate (h1a, s1a)
      end do
      deallocate (eip1)
      if (dh) then
         ! `make_h1_atom` differentiates only integrals; the Kohn-Sham mean
         ! field also keeps part of itself on the quadrature grid, and the
         ! coupled-perturbed right-hand side needs that derivative too.
         call xc_potential_deriv(xc, mol, density, h1, error)
         if (error%has_error()) return
         call solve_mo1_batch(mol, coeff, orbital_energies, n_occ, h1, s1, mo1, &
                              error, xc=xc, reference=density, k_scale=kf, &
                              max_iter=200, tol=use_tol)
      else
         call solve_mo1_batch(mol, coeff, orbital_energies, n_occ, h1, s1, mo1, &
                              error, max_iter=200, tol=use_tol)
      end if
      if (error%has_error()) return
      deallocate (h1, s1)

      ! ---- group 2: the orbital response ------------------------------------
      if (n_frozen > 0) then
         call mp2_orbital_response_term(mo1, sx, fx, xov_st, i2_st, orb, pf=pf_st)
         deallocate (pf_st)
      else
         call mp2_orbital_response_term(mo1, sx, fx, xov_st, i2_st, orb)
      end if
      deallocate (xov_st, i2_st)

      ! ---- group 3: the 2n+1 density response -------------------------------
      call mp2_perturbed_response(mol, coeff, orbital_energies, n_occ, n_frozen, &
                                  fx, sx, erix, mo1, t2, eri_mo, l_mo, gam, &
                                  dm1mo, dt2, ddrel, di, error, tol=use_tol, &
                                  xc=xc, scf_density=density)
      if (error%has_error()) return

      ! The perturbative term's dipole derivative, here because this is the last
      ! point at which `mo1`, `sx`, `fx` and `l_mo` are all still alive -- the
      ! orbital rotations need them, and the line below frees three of the four.
      if (present(dipole_derivatives)) then
         if (.not. dh) then
            call error%set(ERROR_VALIDATION, "the perturbative term's dipole "// &
                           "derivative was asked for without a functional. Over a "// &
                           "Hartree-Fock reference an MP2 dipole derivative is a "// &
                           "different quantity, and this returns the double "// &
                           "hybrid's or nothing.")
            return
         end if
         call correlation_dipole_derivatives(mol, coeff, orbital_energies, n_occ, &
                                             n_frozen, dm1mo, ddrel, mo1, sx, fx, &
                                             l_mo, cscale, dipole_derivatives, error)
         if (error%has_error()) return
      end if

      deallocate (mo1, eri_mo, l_mo, gam, imat)

      ! ---- fold pass 2 into the correlation block ---------------------------
      !
      ! pycc's `_resp(ix, iy)`: the density responses carry the second
      ! perturbation `Y` and the first-order skeletons the first `X`. The
      ! cumulant's response is `mp2_cumulant_2pdm` at `d_Y t2`, built one
      ! perturbation at a time so no `nmo^4 x 3N` stack ever exists.
      do iy = 1, n_pert
         by = (iy - 1)/3 + 1
         cy = iy - 3*(by - 1)
         call mp2_cumulant_2pdm(dt2(:, :, :, :, iy), n_frozen, n_occ, n_mo, dgam)
         do ix = 1, n_pert
            ax = (ix - 1)/3 + 1
            cx = ix - 3*(ax - 1)
            resp = sum(ddrel(:, :, iy)*fx(:, :, ix)) &
                   + sum(dgam*erix(:, :, :, :, ix)) &
                   + sum(di(:, :, iy)*sx(:, :, ix)) &
                   + orb(ix, iy)
            hess_corr(cx, cy, ax, by) = hess_corr(cx, cy, ax, by) + resp
         end do
         deallocate (dgam)
      end do

      ! Once, here, and nowhere else -- the gradient's own rule, one order up.
      ! `hess_ref` is deliberately not scaled: it is the reference's skeleton,
      ! and a double-hybrid caller discards it anyway (header).
      if (dh) hess_corr = cscale*hess_corr

      deallocate (gradient, dm1mo, w_ao, t2, fx, sx, erix, orb, dt2, ddrel, di)
      if (allocated(vk_drel)) deallocate (vk_drel, kd_drel_mo)
      if (allocated(vkpf)) deallocate (vkpf)
   end subroutine mp2_correlation_hessian

   subroutine mp2_ks_first_order_fock(mol, xc, scf_density, coeff, n_occ, erix, &
                                      k_scale, fx, error)
      !! Lift the Hartree-Fock first-order skeleton Fock to the Kohn-Sham one
      !!
      !!     f^(X)_pq -> f^(X)_pq + (1 - k) sum_m <pm|mq>^(X) + (C^T V_xc^(X) C)_pq
      !!
      !! `mp2_first_order_skeletons` builds the Hartree-Fock skeleton and its
      !! text is part of the no-`xc` path's bit-for-bit contract (module
      !! header), so the lift lives here instead: the exchange trace is
      !! recovered from the `erix` stack that builder already returns --
      !! `<pm|mq>^(X)` is its fold's `chem(p, m, m, q)` said in physicist
      !! order -- and the potential's skeleton derivative comes from one
      !! `xc_potential_deriv` grid pass covering every perturbation, matrices
      !! and grid held fixed. `test_mqc_mp2_hessian_skeleton_fock` differences
      !! the lifted matrix against the Kohn-Sham Fock it derives.
      type(czt_molecule_t), intent(in) :: mol
      type(xc_context_t), intent(inout) :: xc
      real(dp), intent(in) :: scf_density(:, :)   !! Converged reference density
      real(dp), intent(in) :: coeff(:, :)
      integer, intent(in) :: n_occ                !! Full doubly-occupied count
      real(dp), intent(in) :: erix(:, :, :, :, :)  !! `<pq|rs>^(X)`, physicist
      real(dp), intent(in) :: k_scale             !! The functional's fraction
      real(dp), intent(inout) :: fx(:, :, :)      !! Lifted in place
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: vxc1(:, :, :, :)
      real(dp) :: acc, cx
      integer :: n_ao, n_mo, natm, a, comp, x, p, q, m

      if (error%has_error()) return
      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)
      natm = mol%natm
      allocate (vxc1(n_ao, n_ao, 3, natm))
      vxc1 = 0.0_dp
      call xc_potential_deriv(xc, mol, scf_density, vxc1, error)
      if (error%has_error()) return

      cx = 1.0_dp - k_scale
      do a = 1, natm
         do comp = 1, 3
            x = 3*(a - 1) + comp
            do q = 1, n_mo
               do p = 1, n_mo
                  acc = 0.0_dp
                  do m = 1, n_occ
                     acc = acc + erix(p, m, m, q, x)
                  end do
                  fx(p, q, x) = fx(p, q, x) + cx*acc
               end do
            end do
            fx(:, :, x) = fx(:, :, x) &
                          + matmul(transpose(coeff), &
                                   matmul(vxc1(:, :, comp, a), coeff))
         end do
      end do
      deallocate (vxc1)
   end subroutine mp2_ks_first_order_fock

   subroutine ks_occ_fold_patch(t, d, n_occ, k_scale, vk, mat)
      !! Lift a Lagrangian's occupied-column fold to the Kohn-Sham reference
      !!
      !! `mp2_mo_lagrangian` and `mp2_skeleton_lagrangian` fold
      !! `D_rs (L_rpsq + L_rqsp)` at full exchange and close with `-1/2`;
      !! over a Kohn-Sham reference the same fold carries the functional's
      !! exchange fraction and the kernel. Both fixes are additive on top of
      !! the Hartree-Fock result,
      !!
      !!     mat(p, q) += -1/2 [ (1-k) sum_rs D(r,s) (T(r,p,q,s) + T(r,q,p,s))
      !!                         + vk(p, q) ],   q occupied
      !!
      !! and living here rather than inside those routines keeps their text
      !! -- and therefore the no-`xc` path's bits -- untouched (module
      !! header). The rescale stays a threaded loop for the folds' own
      !! reason: `p` and `q` are strided in `T`, and it is one power of
      !! `n_mo` cheaper than the gemm terms beside it.
      real(dp), intent(in) :: t(:, :, :, :)   !! `<pq|rs>` or `<pq|rs>^(X)`
      real(dp), intent(in) :: d(:, :)         !! The fold's density
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: k_scale
      real(dp), intent(in) :: vk(:, :)
         !! The kernel share: `ks_pair_kernel` at `d` for the unperturbed
         !! fold, the matching `ks_kernel_deriv_stack` slice for a skeleton-
         !! perturbed one.
      real(dp), intent(inout) :: mat(:, :)

      real(dp) :: acc, cx
      integer :: n_mo, p, q, r, s

      n_mo = size(d, 1)
      cx = 1.0_dp - k_scale
      !$omp parallel do collapse(2) default(none) schedule(static) &
      !$omp    shared(mat, d, t, vk, n_mo, n_occ, cx) private(p, q, r, s, acc)
      do q = 1, n_occ
         do p = 1, n_mo
            acc = 0.0_dp
            do s = 1, n_mo
               do r = 1, n_mo
                  acc = acc + d(r, s)*(t(r, p, q, s) + t(r, q, p, s))
               end do
            end do
            mat(p, q) = mat(p, q) - 0.5_dp*(cx*acc + vk(p, q))
         end do
      end do
      !$omp end parallel do
   end subroutine ks_occ_fold_patch

   subroutine ks_pair_kernel(xc, mol, coeff, scf_density, m_mo, vout, error)
      !! The exchange-correlation kernel's share of one pair contraction
      !!
      !! Every kernel site on this ladder has the shape
      !! `test_mqc_mp2_hessian_ks_operator` pins,
      !!
      !!     sum_rs M(r,s) (L(r,p,s,q) + L(r,q,s,p))
      !!
      !! and its kernel share is the reference operator applied to a symmetric
      !! generalized density -- never a four-index `f_xc` tensor:
      !!
      !!     vout = 4 C^T f_xc[ C M_sym C^T ] C,   M_sym = (M + M^T)/2
      !!
      !! The factor is anchored to the coupled-perturbed operator this must
      !! agree with: `response_operator` applies `2 G[Dt]` at
      !! `Dt = 2 C M_sym C^T`, and `xc_kernel_apply` is linear in its trial
      !! density. The contraction annihilates the antisymmetric part of `M`
      !! -- the operator test's deliberately asymmetric trial -- so callers
      !! hand over one-sided densities (a vo rotation, one triangle of a
      !! core<->active divide) exactly as they hold them.
      type(xc_context_t), intent(inout) :: xc
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: scf_density(:, :)   !! Where the kernel lives
      real(dp), intent(in) :: m_mo(:, :)          !! Generalized density, MO
      real(dp), intent(out) :: vout(:, :)         !! (n_mo, n_mo)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: msym(:, :), m_ao(:, :), vao(:, :)
      integer :: n_ao, n_mo

      if (error%has_error()) return
      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)
      allocate (msym(n_mo, n_mo), m_ao(n_ao, n_ao), vao(n_ao, n_ao))
      msym = 0.5_dp*(m_mo + transpose(m_mo))
      m_ao = matmul(coeff, matmul(msym, transpose(coeff)))
      vao = 0.0_dp
      call xc_kernel_apply(xc, mol, scf_density, m_ao, vao, error)
      if (error%has_error()) return
      vout = 4.0_dp*matmul(transpose(coeff), matmul(vao, coeff))
      deallocate (msym, m_ao, vao)
   end subroutine ks_pair_kernel

   subroutine ks_kernel_deriv_stack(xc, mol, coeff, scf_density, m_mo, kd_mo, error)
      !! The skeleton nuclear derivative of `ks_pair_kernel`, every perturbation
      !!
      !!     kd_mo(:, :, x) = 4 C^T ( d/dX f_xc-apply[C M_sym C^T] ) C
      !!
      !! at fixed matrices and fixed grid -- `xc_kernel_deriv`, which already
      !! carries both families of skeleton motion (the basis functions under
      !! the integral and the densities on the grid, `g_xc` included), so one
      !! grid pass covers all `3 natm` slices. What it deliberately does not
      !! carry -- the MO rotations and the reference density's matrix
      !! response -- is `ks_dl_pair_kernel`'s per-perturbation business.
      type(xc_context_t), intent(inout) :: xc
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: scf_density(:, :)
      real(dp), intent(in) :: m_mo(:, :)
      real(dp), allocatable, intent(out) :: kd_mo(:, :, :)   !! (n_mo, n_mo, 3*natm)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: msym(:, :), m_ao(:, :), h1(:, :, :, :)
      integer :: n_ao, n_mo, natm, a, comp, x

      if (error%has_error()) return
      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)
      natm = mol%natm
      allocate (msym(n_mo, n_mo), m_ao(n_ao, n_ao))
      msym = 0.5_dp*(m_mo + transpose(m_mo))
      m_ao = matmul(coeff, matmul(msym, transpose(coeff)))
      allocate (h1(n_ao, n_ao, 3, natm))
      h1 = 0.0_dp
      call xc_kernel_deriv(xc, mol, scf_density, m_ao, h1, error)
      if (error%has_error()) return
      allocate (kd_mo(n_mo, n_mo, 3*natm))
      do a = 1, natm
         do comp = 1, 3
            x = 3*(a - 1) + comp
            kd_mo(:, :, x) = 4.0_dp*matmul(transpose(coeff), &
                                           matmul(h1(:, :, comp, a), coeff))
         end do
      end do
      deallocate (msym, m_ao, h1)
   end subroutine ks_kernel_deriv_stack

   subroutine ks_dl_pair_kernel(xc, mol, coeff, scf_density, m_mo, vk_m, kd_x, &
                                u, ddmat_ao, w, error)
      !! The kernel share of a full-derivative fold: `d_X L` against `M`
      !!
      !! What `ks_pair_kernel` is to `L`, this is to the derivative weight
      !! `dl` -- the total `d_X` of `4 C^T f_xc[C M_sym C^T] C` at fixed `M`,
      !! in the four pieces the product rule gives:
      !!
      !!     w = U^T Vk + Vk U                    the free MO pair rotating
      !!       + kd_x                             the grid skeleton (stack slice)
      !!       + apply[U M_sym + M_sym U^T]       the contracted density rotating
      !!       + 4 C^T g_xc[dD, C M_sym C^T] C    the kernel's coefficients
      !!                                          moving with the reference
      !!                                          density's matrix response dD
      !!
      !! The `g_xc` term is `xc_kernel2_apply`; the skeleton share of the
      !! density's motion is already inside `kd_x`, so `ddmat_ao` is the
      !! matrix response alone -- `nuclear_apply`'s first-order density,
      !! `2 (C U_occ C_occ^T + h.c.)`.
      type(xc_context_t), intent(inout) :: xc
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: scf_density(:, :)
      real(dp), intent(in) :: m_mo(:, :)      !! The fixed generalized density
      real(dp), intent(in) :: vk_m(:, :)      !! `ks_pair_kernel` at `m_mo`
      real(dp), intent(in) :: kd_x(:, :)      !! One slice of the deriv stack
      real(dp), intent(in) :: u(:, :)         !! `U^X`, full MO rotation
      real(dp), intent(in) :: ddmat_ao(:, :)  !! Reference density's response
      real(dp), intent(out) :: w(:, :)        !! (n_mo, n_mo)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: msym(:, :), rot(:, :), tmp(:, :)
      real(dp), allocatable :: m_ao(:, :), vao(:, :)
      integer :: n_ao, n_mo

      if (error%has_error()) return
      n_ao = size(coeff, 1)
      n_mo = size(coeff, 2)
      allocate (msym(n_mo, n_mo), rot(n_mo, n_mo), tmp(n_mo, n_mo))
      msym = 0.5_dp*(m_mo + transpose(m_mo))

      w = matmul(transpose(u), vk_m) + matmul(vk_m, u) + kd_x

      rot = matmul(u, msym) + matmul(msym, transpose(u))
      call ks_pair_kernel(xc, mol, coeff, scf_density, rot, tmp, error)
      if (error%has_error()) return
      w = w + tmp

      allocate (m_ao(n_ao, n_ao), vao(n_ao, n_ao))
      m_ao = matmul(coeff, matmul(msym, transpose(coeff)))
      vao = 0.0_dp
      call xc_kernel2_apply(xc, mol, scf_density, ddmat_ao, m_ao, vao, error)
      if (error%has_error()) return
      w = w + 4.0_dp*matmul(transpose(coeff), matmul(vao, coeff))
      deallocate (msym, rot, tmp, m_ao, vao)
   end subroutine ks_dl_pair_kernel

   subroutine correlation_dipole_derivatives(mol, coeff, orbital_energies, n_occ, &
                                             n_frozen, dm1mo, ddrel, mo1, sx, fx, &
                                             l_mo, scale, ddip_dr, error)
      !! `c d mu_corr / dR`, from the relaxed density and its response
      !!
      !! The correlation term reaches the dipole through one object, the relaxed
      !! density, so its derivative has two terms and no others:
      !!
      !!     d mu_a / dR_(X,b) = - Tr(dD_rel/dR_(X,b)  r_a)
      !!                         - Tr(D_rel  d r_a / dR_(X,b))
      !!
      !! There is no nuclear term: the nuclei belong to the reference's dipole,
      !! the same split the perturbative gradient makes when it leaves nuclear
      !! repulsion to `czt_scf_gradient`.
      !!
      !! **`ddrel` is not the AO derivative**, and assuming it was is the
      !! mistake this routine exists to not make. It is the derivative of the
      !! relaxed density's *matrix elements*, while the dipole contracts the AO
      !! density `C D C^T`, whose derivative also moves the orbitals:
      !!
      !!     dD_AO/dR = C [ dD_MO/dR + U D_MO + D_MO U^T ] C^T
      !!
      !! The first reading of this omitted the two rotation terms and was wrong
      !! by 1.11e-02 -- the same order as the correlation dipole itself, 1.26e-02
      !! -- while remaining smooth and obeying the translational sum rule.
      !! `test_mqc_dh_dipole_deriv` differences the dipole and is what caught it.
      !!
      !! `U` is the full `nmo x nmo` rotation rather than the occupied columns
      !! the response solves for, because the relaxed density has virtual-virtual
      !! and occupied-occupied blocks that those columns never touch --
      !! `mp2_full_u` assembles it, with the frozen-core branch it documents.
      type(czt_molecule_t), intent(in), target :: mol
      real(dp), intent(in) :: coeff(:, :)
      real(dp), intent(in) :: orbital_energies(:)
      integer, intent(in) :: n_occ, n_frozen
      real(dp), intent(in) :: dm1mo(:, :)        !! Relaxed correlation density, MO
      real(dp), intent(in) :: ddrel(:, :, :)     !! Its matrix elements' derivative
      real(dp), intent(in) :: mo1(:, :, :, :)    !! (n_mo, n_occ, 3, natm)
      real(dp), intent(in) :: sx(:, :, :), fx(:, :, :)   !! S^(X), f^(X), MO
      real(dp), intent(in) :: l_mo(:, :, :, :)
      real(dp), intent(in) :: scale              !! The functional's PT2 coefficient
      real(dp), allocatable, intent(out) :: ddip_dr(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dip(:, :, :), ddip(:, :, :, :, :), drel_ao(:, :)
      real(dp), allocatable :: dd_ao(:, :), dd_mo(:, :), u(:, :)
      real(dp) :: origin(3)
      integer :: natm, ia, a, b, y, n_ao

      natm = mol%natm
      n_ao = size(coeff, 1)

      ! The nuclear centroid, matching `assemble_dipole_derivatives`: the caller
      ! adds the two halves, so they have to be expanded about the same point
      ! even though the total does not depend on it.
      origin = 0.0_dp
      do ia = 1, natm
         origin = origin + mol%coords(:, ia)
      end do
      if (natm > 0) origin = origin/real(natm, dp)

      call multipole_matrices(mol, origin, 1, dip, error)
      if (error%has_error()) return
      call dipole_integral_derivatives(mol, origin, ddip, error)
      if (error%has_error()) return

      drel_ao = matmul(coeff, matmul(dm1mo, transpose(coeff)))
      allocate (ddip_dr(3, 3*natm), dd_ao(n_ao, n_ao))
      allocate (dd_mo(size(dm1mo, 1), size(dm1mo, 2)))
      ddip_dr = 0.0_dp

      do ia = 1, natm
         do b = 1, 3
            y = 3*(ia - 1) + b
            call mp2_full_u(mo1(:, :, b, ia), sx(:, :, y), n_occ, u, &
                            n_frozen=n_frozen, fx1=fx(:, :, y), l_mo=l_mo, &
                            orbital_energies=orbital_energies)
            dd_mo = ddrel(:, :, y) + matmul(u, dm1mo) + matmul(dm1mo, transpose(u))
            dd_ao = matmul(coeff, matmul(dd_mo, transpose(coeff)))
            do a = 1, 3
               ddip_dr(a, y) = -scale*(sum(dd_ao*dip(:, :, a)) &
                                       + sum(drel_ao*ddip(:, :, a, b, ia)))
            end do
            deallocate (u)
         end do
      end do

      deallocate (dip, ddip, drel_ao, dd_ao, dd_mo)
   end subroutine correlation_dipole_derivatives

end module mqc_czt_mp2_hessian
