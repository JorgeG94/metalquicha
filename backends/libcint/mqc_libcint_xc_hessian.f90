!! The exchange-correlation contribution to an analytic nuclear Hessian
module mqc_libcint_xc_hessian
   !! d2 E_xc / dR dR at fixed molecular orbital coefficients
   !!
   !! The explicit half of a Kohn-Sham Hessian: what the exchange-correlation
   !! energy does when the nuclei move and the density matrix is held still.
   !! The other half -- the density's own response to the perturbation -- is the
   !! coupled-perturbed problem and belongs with the rest of the response terms.
   !!
   !! **The grid does not respond here, deliberately.** The quadrature is
   !! atom-centred, so a rigorous second derivative would carry the points
   !! moving with their owner and the partition weights being reweighted by
   !! every nucleus. Both are omitted, which is what PySCF's `Hessian.grid_response`
   !! defaults to as well -- so this is comparable against it term for term
   !! rather than approximately. The omission is a real approximation and it is
   !! not free: the same term left out of a *gradient* costs about 1e-4, which
   !! is why `xc_gradient` next door includes it. It is left out here because
   !! the second derivative of a Becke weight is a `(3, natm, 3, natm)` object
   !! per grid point that cannot be stored and is expensive to contract, and
   !! because the reference everyone checks against omits it too. See
   !! `becke_cutoff_second_derivative`, which is the groundwork for adding it.
   !!
   !! **What the two terms are.** With `rho = sum_uv D_uv chi_u chi_v` and a
   !! basis function's nuclear derivative being minus its position derivative,
   !!
   !!     d rho / dA_c   = -2 sum_{u in A, v} D_uv (d_c chi_u) chi_v
   !!     d2 rho / dA_c dB_d
   !!         = 2 delta_AB sum_{u in A, v} D_uv (d_c d_d chi_u) chi_v
   !!         + 2 sum_{u in A, v in B} D_uv (d_c chi_u) (d_d chi_v)
   !!
   !! and the energy's second derivative is the kernel against the first
   !! derivatives plus the potential against the second:
   !!
   !!     d2 E / dA dB = sum_g w_g [ f_rr (d rho/dA)(d rho/dB)
   !!                                + v_rho (d2 rho/dA dB) ]
   !!
   !! For a GGA every one of those gains a `sigma` channel, which brings a third
   !! derivative of the basis functions: `sigma` already costs one position
   !! derivative to form, and two nuclear derivatives on top of it land on
   !! `d_c d_d d_e chi`. That is why this refuses a GGA rather than quietly
   !! omitting the terms -- `eval_ao_block` produces two derivatives and the
   !! third is a separate piece of work.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_ao, only: eval_ao_block, shell_extents, block_significant_aos, &
                             AO_DERIV3_COMP, eval_rho
   use mqc_libcint_xc, only: xc_context_t, xc_grid_kernel_quantities, &
                             ensure_nlc_grid
   use mqc_libcint_vv10, only: vv10_nlc, vv10_hessian_kernel
   implicit none
   private

   public :: xc_hessian
   public :: xc_potential_hessian
   public :: xc_gradient_fixed_grid
   public :: xc_potential_deriv
   public :: vv10_hessian
   public :: vv10_potential_deriv

contains

   subroutine xc_hessian(ctx, mol, density, hess, error)
      !! Accumulate the exchange-correlation second derivatives into `hess`
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: hess(:, :, :, :)   !! (3, 3, natm, natm), in place
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), rho_grad(:, :), vrho(:), vsigma(:)
      real(dp), allocatable :: frr(:), frs(:), fss(:)
      real(dp), allocatable :: tau(:), vtau(:), frt(:), fst(:), ftt(:)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), dchi(:, :)
      real(dp), allocatable :: dmu(:, :), zmat(:, :), ymat(:, :), wg(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:)
      integer, allocatable :: c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig
      integer :: ia, ja, a, b, i, j, ig, p, q, c, d, e, ih
      real(dp) :: acc
      logical :: gga, mgga
      real(dp), allocatable :: ao_d3(:, :, :), dmu_d(:, :, :)
      real(dp), allocatable :: dgrad(:, :, :), dsig(:, :), wsig(:, :), dtau(:, :)
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
      integer, parameter :: TRIP(3, 3, 3) = reshape( &
                            [1, 2, 3, 2, 4, 5, 3, 5, 6, &
                             2, 4, 5, 4, 7, 8, 5, 8, 9, &
                             3, 5, 6, 5, 8, 9, 6, 9, 10], [3, 3, 3])
         !! Which packed third-derivative component a Cartesian triple is.

      if (error%has_error()) return

      mgga = ctx%any_mgga
      gga = ctx%any_gga .or. mgga

      natm = size(mol%coords, 2)
      nao = mol%nao
      npts = size(ctx%grid%weights)

      ! rho, the potential and the kernel over the whole grid. The GGA channels
      ! come back too and are unused here, which the refusal above makes safe.
      if (mgga) then
         call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                        frr, frs, fss, error, tau=tau, vtau=vtau, &
                                        frt=frt, fst=fst, ftt=ftt)
      else
         call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                        frr, frs, fss, error)
      end if
      if (error%has_error()) return

      allocate (extents(mol%nbas))
      call shell_extents(mol, ctx%screen_tol, extents)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (c_offsets(natm), c_counts(natm))
      allocate (d_sig(nao, nao))

      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)
         if (n_sig == 0) cycle

         do jsig = 1, n_sig
            do isig = 1, n_sig
               d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
            end do
         end do

         if (gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                               grad=ao_grad, hess=ao_hess, deriv3=ao_d3, &
                               shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                               grad=ao_grad, hess=ao_hess, shell_mask=shell_mask, &
                               ao_offset=ao_offset, n_ao_out=n_sig)
         end if
         if (error%has_error()) return

         if (allocated(wg)) deallocate (wg)
         allocate (wg(nb))
         wg = ctx%grid%weights(g0:g1)

         ! (D chi)(u, g), one gemm, shared by every term below.
         if (allocated(dmu)) deallocate (dmu)
         allocate (dmu(nb, n_sig))
         call pic_gemm(ao(1:nb, 1:n_sig), d_sig(1:n_sig, 1:n_sig), dmu, beta=0.0_dp)

         ! d rho / dA_c at every point of the block, as (3*natm, nb).
         if (allocated(dchi)) deallocate (dchi)
         allocate (dchi(3*natm, nb))
         dchi = 0.0_dp
         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               do ig = 1, nb
                  acc = 0.0_dp
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     acc = acc + ao_grad(ig, i, a)*dmu(ig, i)
                  end do
                  dchi(3*(ia - 1) + a, ig) = -2.0_dp*acc
               end do
            end do
         end do

         ! d(grad rho)_d / dA_c, and from it d sigma / dA_c. Both are needed by
         ! the kernel terms and by the potential term below, so they are formed
         ! once here rather than twice.
         if (gga) then
            if (allocated(dmu_d)) deallocate (dmu_d)
            allocate (dmu_d(nb, n_sig, 3))
            do d = 1, 3
               call pic_gemm(ao_grad(1:nb, 1:n_sig, d), d_sig(1:n_sig, 1:n_sig), &
                             dmu_d(:, :, d), beta=0.0_dp)
            end do

            if (allocated(dgrad)) deallocate (dgrad)
            if (allocated(dsig)) deallocate (dsig)
            allocate (dgrad(3*natm, 3, nb), dsig(3*natm, nb))
            dgrad = 0.0_dp
            dsig = 0.0_dp
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               do c = 1, 3
                  do d = 1, 3
                     do ig = 1, nb
                        acc = 0.0_dp
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + ao_hess(ig, i, PAIR(c, d))*dmu(ig, i) &
                                 + ao_grad(ig, i, c)*dmu_d(ig, i, d)
                        end do
                        dgrad(3*(ia - 1) + c, d, ig) = -2.0_dp*acc
                     end do
                  end do
               end do
            end do
            do ig = 1, nb
               do p = 1, 3*natm
                  acc = 0.0_dp
                  do d = 1, 3
                     acc = acc + rho_grad(g0 + ig - 1, d)*dgrad(p, d, ig)
                  end do
                  dsig(p, ig) = 2.0_dp*acc
               end do
            end do
         end if

         ! d tau / dA. With tau = 1/2 sum_d sum_uv D_uv (d_d chi_u)(d_d chi_v),
         ! one nuclear derivative lands on a second position derivative and the
         ! half cancels against the two equal halves of the product rule.
         if (mgga) then
            if (allocated(dtau)) deallocate (dtau)
            allocate (dtau(3*natm, nb))
            dtau = 0.0_dp
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               do c = 1, 3
                  do ig = 1, nb
                     acc = 0.0_dp
                     do d = 1, 3
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + ao_hess(ig, i, PAIR(c, d))*dmu_d(ig, i, d)
                        end do
                     end do
                     dtau(3*(ia - 1) + c, ig) = -acc
                  end do
               end do
            end do
         end if

         ! Kernel terms: the second functional derivatives against the first
         ! derivatives of their own arguments. For an LDA that is `f_rr` alone;
         ! a GGA adds the two cross terms and `f_ss`.
         do ig = 1, nb
            do p = 1, 3*natm
               do q = 1, 3*natm
                  acc = frr(g0 + ig - 1)*dchi(p, ig)*dchi(q, ig)
                  if (gga) then
                     acc = acc + frs(g0 + ig - 1) &
                           *(dchi(p, ig)*dsig(q, ig) + dsig(p, ig)*dchi(q, ig)) &
                           + fss(g0 + ig - 1)*dsig(p, ig)*dsig(q, ig)
                  end if
                  if (mgga) then
                     acc = acc + frt(g0 + ig - 1) &
                           *(dchi(p, ig)*dtau(q, ig) + dtau(p, ig)*dchi(q, ig)) &
                           + fst(g0 + ig - 1) &
                           *(dsig(p, ig)*dtau(q, ig) + dtau(p, ig)*dsig(q, ig)) &
                           + ftt(g0 + ig - 1)*dtau(p, ig)*dtau(q, ig)
                  end if
                  hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) = &
                     hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) &
                     + wg(ig)*acc
               end do
            end do
         end do

         ! Potential term against the second derivative of sigma, first half:
         ! the product of two first derivatives of grad rho. Same shape as the
         ! kernel terms, so it goes with them.
         if (gga) then
            do ig = 1, nb
               do p = 1, 3*natm
                  do q = 1, 3*natm
                     acc = 0.0_dp
                     do d = 1, 3
                        acc = acc + dgrad(p, d, ig)*dgrad(q, d, ig)
                     end do
                     hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) = &
                        hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) &
                        + 2.0_dp*wg(ig)*vsigma(g0 + ig - 1)*acc
                  end do
               end do
            end do
         end if

         ! Potential term, second piece: the off-diagonal product of two first
         ! derivatives of different functions. Built as nine weighted overlaps
         ! `Z_cd(u,v) = sum_g w v_rho (d_c chi_u)(d_d chi_v)` and then masked
         ! onto atom pairs, which keeps the grid sum a gemm and the atom
         ! bookkeeping out of the inner loop.
         if (allocated(zmat)) deallocate (zmat)
         allocate (zmat(n_sig, n_sig))
         do a = 1, 3
            do b = 1, 3
               call weighted_overlap(ao_grad(1:nb, 1:n_sig, a), ao_grad(1:nb, 1:n_sig, b), &
                                     wg*vrho(g0:g1), zmat)
               do ja = 1, natm
                  if (c_counts(ja) == 0) cycle
                  do ia = 1, natm
                     if (c_counts(ia) == 0) cycle
                     acc = 0.0_dp
                     do j = c_offsets(ja) + 1, c_offsets(ja) + c_counts(ja)
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + d_sig(i, j)*zmat(i, j)
                        end do
                     end do
                     hess(a, b, ia, ja) = hess(a, b, ia, ja) + 2.0_dp*acc
                  end do
               end do
            end do
         end do

         ! Potential term against the second derivative of sigma, second half:
         ! `grad rho` contracted with `d2(grad rho)`. This is the only place the
         ! third derivatives of the basis functions are used, and they appear
         ! only on the diagonal atom block -- both nuclear derivatives landing
         ! on the same function is what produces a third position derivative.
         if (gga) then
            if (allocated(wsig)) deallocate (wsig)
            allocate (wsig(nb, 3))
            do d = 1, 3
               wsig(:, d) = 2.0_dp*wg*vsigma(g0:g1)*rho_grad(g0:g1, d)
            end do
            if (allocated(zmat)) deallocate (zmat)
            allocate (zmat(n_sig, n_sig))

            do d = 1, 3
               do c = 1, 3
                  do e = 1, 3
                     ! Diagonal: d3 chi against chi, and d2 chi against d chi.
                     call weighted_overlap(ao_d3(1:nb, 1:n_sig, TRIP(c, e, d)), &
                                           ao(1:nb, 1:n_sig), wsig(:, d), zmat)
                     call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                     natm, n_sig, c, e, 2.0_dp, diagonal=.true.)
                     call weighted_overlap(ao_hess(1:nb, 1:n_sig, PAIR(c, e)), &
                                           ao_grad(1:nb, 1:n_sig, d), wsig(:, d), zmat)
                     call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                     natm, n_sig, c, e, 2.0_dp, diagonal=.true.)
                     ! Off-diagonal: one derivative on each centre.
                     call weighted_overlap(ao_hess(1:nb, 1:n_sig, PAIR(c, d)), &
                                           ao_grad(1:nb, 1:n_sig, e), wsig(:, d), zmat)
                     call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                     natm, n_sig, c, e, 2.0_dp, diagonal=.false.)
                     call weighted_overlap(ao_grad(1:nb, 1:n_sig, c), &
                                           ao_hess(1:nb, 1:n_sig, PAIR(e, d)), wsig(:, d), zmat)
                     call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                     natm, n_sig, c, e, 2.0_dp, diagonal=.false.)
                  end do
               end do
            end do
         end if

         ! `v_tau` against the second derivative of tau. Same shape as the sigma
         ! term above: a diagonal piece where both nuclear derivatives land on
         ! one centre and reach a third position derivative, and an off-diagonal
         ! piece with one on each.
         if (mgga) then
            if (allocated(wsig)) deallocate (wsig)
            allocate (wsig(nb, 1))
            wsig(:, 1) = wg*vtau(g0:g1)
            if (allocated(zmat)) deallocate (zmat)
            allocate (zmat(n_sig, n_sig))
            do c = 1, 3
               do e = 1, 3
                  do d = 1, 3
                     call weighted_overlap(ao_d3(1:nb, 1:n_sig, TRIP(c, e, d)), &
                                           ao_grad(1:nb, 1:n_sig, d), wsig(:, 1), zmat)
                     call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                     natm, n_sig, c, e, 1.0_dp, diagonal=.true.)
                     call weighted_overlap(ao_hess(1:nb, 1:n_sig, PAIR(c, d)), &
                                           ao_hess(1:nb, 1:n_sig, PAIR(e, d)), &
                                           wsig(:, 1), zmat)
                     call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                     natm, n_sig, c, e, 1.0_dp, diagonal=.false.)
                  end do
               end do
            end do
         end if

         ! Potential term, first piece: both derivatives on the same function,
         ! so it lands on the diagonal atom block only.
         if (allocated(ymat)) deallocate (ymat)
         allocate (ymat(n_sig, n_sig))
         do p = 1, 6
            call hess_component(p, a, b)
            call weighted_overlap(ao_hess(1:nb, 1:n_sig, p), ao(1:nb, 1:n_sig), &
                                  wg*vrho(g0:g1), ymat)
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               acc = 0.0_dp
               do j = 1, n_sig
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     acc = acc + d_sig(i, j)*ymat(i, j)
                  end do
               end do
               hess(a, b, ia, ia) = hess(a, b, ia, ia) + 2.0_dp*acc
               if (a /= b) hess(b, a, ia, ia) = hess(b, a, ia, ia) + 2.0_dp*acc
            end do
         end do
      end do
   end subroutine xc_hessian

   subroutine xc_potential_hessian(ctx, mol, density, pmat, hess, error)
      !! `d2/dRdR Tr(P V_xc[D])` with both densities and the grid held fixed,
      !! accumulated in place
      !!
      !! The second derivative of the linear form `xc_potential_gradient`
      !! differentiates once: the pass-one term of a double-hybrid Hessian
      !! contracts the relaxed density against `V_xc^(XY)`, because the
      !! reference operator's second derivative contains the potential's. As
      !! there, `P` is symmetric, indefinite, integrates to zero, and appears
      !! only linearly -- the functional is evaluated at `D` alone.
      !!
      !! Writing `u = (rho, sigma)` for the reference's variables and `u_P` for
      !! the linearisation in `P` -- `rho_P`, and
      !! `sigma_P = 2 grad rho . grad rho_P` -- the integrand collects into
      !! four groups by how far up the functional-derivative ladder each sits:
      !!
      !!     g_xc u^(X) u^(Y) u_P                       the kernel moving twice
      !!     f_xc (u^(X) u_P^(Y) + u_P^(X) u^(Y))       once, against P's motion
      !!     f_xc u_P u^(XY)                            once, against d2(rho, sigma)
      !!     v_xc u_P^(XY)                              the potential, against
      !!                                                d2(rho_P, sigma_P)
      !!
      !! every product summed over both channels -- the expansion
      !! `xc_kernel_apply` documents for `f_xc`, one order up. All densities
      !! are AO motion only: the matrices never respond, and neither does the
      !! grid, so the check differences the fixed-grid first derivative of the
      !! same linear form -- `xc_potential_deriv` contracted with `P`, the
      !! same approximation on both sides. Not `xc_potential_gradient` with
      !! `fixed_grid`, though that flag exists for this comparison: its
      !! deposits never forward `moving=.false.` to `accumulate_channel`, so
      !! that branch keeps the grid-point-motion term while dropping the
      !! partition weights, and the result is the exact derivative of nothing
      !! -- the test's helper carries the measurement.
      !!
      !! **A separate routine, deliberately, where the double-hybrid gradient
      !! argued itself into a merge.** That merge (`libcint_mp2_gradient`) was
      !! between two calculations that differ in four places and agree
      !! everywhere else -- one body, four branches. Here there is no argument
      !! at which this coincides with `xc_hessian`: handing it `pmat = density`
      !! yields the second derivative of `int w [v_rho rho + 2 v_sigma sigma]`,
      !! which is nothing `xc_hessian` returns, so an optional second density
      !! there would interleave two disjoint integrands behind one name -- the
      !! near-copy risk relocated, not removed. What genuinely repeats is
      !! factored instead: the weighted overlaps, the atom masking, the packed
      !! derivative indices, and the one block that would otherwise have become
      !! a third spelled-out copy -- the `grad . d2(grad)` contraction -- is
      !! `d2grad_sigma_term`, called three times here with three
      !! weight/density pairings where `xc_hessian` inlines its one.
      !!
      !! The last two groups above are `xc_hessian`'s potential machinery with,
      !! respectively, the kernel contracted with `(rho_P, sigma_P)` standing
      !! in for the potential -- the same stand-in `xc_potential_gradient`
      !! calls `coef_rho`, one derivative down -- and `P` standing in for `D`.
      !!
      !! LDA and GGA. A meta-GGA is refused explicitly: it would need the tau
      !! channels of the third functional derivative, which nothing here
      !! provides.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)   !! The converged reference density
      real(dp), intent(in) :: pmat(:, :)      !! The density `V_xc` is traced against
      real(dp), intent(inout) :: hess(:, :, :, :)   !! (3, 3, natm, natm), in place
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), rho_grad(:, :), vrho(:), vsigma(:)
      real(dp), allocatable :: frr(:), frs(:), fss(:)
      real(dp), allocatable :: grrr(:), grrs(:), grss(:), gsss(:)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: ao_d3(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), p_sig(:, :)
      real(dp), allocatable :: dmu(:, :), pmu(:, :)
      real(dp), allocatable :: dmu_d(:, :, :), pmu_d(:, :, :)
      real(dp), allocatable :: dchi(:, :), pchi(:, :)
      real(dp), allocatable :: dgrad(:, :, :), pgrad(:, :, :)
      real(dp), allocatable :: dsig(:, :), psig(:, :)
      real(dp), allocatable :: rho_p(:), p_grad(:, :), sig_p(:)
      real(dp), allocatable :: vrho_eff(:), vsig_eff(:)
      real(dp), allocatable :: wsig(:, :), zmat(:, :), wg(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:)
      integer, allocatable :: c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig
      integer :: ia, a, b, i, ig, gg, p, q, c, d, ih
      real(dp) :: acc, accp, wa, wb, wc
      logical :: gga
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      if (error%has_error()) return

      if (ctx%any_mgga) then
         call error%set(ERROR_VALIDATION, "the exchange-correlation potential Hessian "// &
                        "is LDA and GGA only: a meta-GGA needs the tau channels of the "// &
                        "third functional derivative, which nothing here provides. "// &
                        "Refused rather than computed with those terms missing.")
         return
      end if
      gga = ctx%any_gga

      natm = size(mol%coords, 2)
      nao = mol%nao
      npts = size(ctx%grid%weights)

      ! The whole ladder at once: potential, kernel, and the four third
      ! derivatives, floored together at KERNEL_RHO_FLOOR so no order crosses
      ! a cutoff the others respect.
      call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                     frr, frs, fss, error, &
                                     grrr=grrr, grrs=grrs, grss=grss, gsss=gsss)
      if (error%has_error()) return

      allocate (extents(mol%nbas))
      call shell_extents(mol, ctx%screen_tol, extents)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (c_offsets(natm), c_counts(natm))
      allocate (d_sig(nao, nao), p_sig(nao, nao))

      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)
         if (n_sig == 0) cycle

         do jsig = 1, n_sig
            do isig = 1, n_sig
               d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
               p_sig(isig, jsig) = pmat(ao_list(isig), ao_list(jsig))
            end do
         end do

         if (gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                               grad=ao_grad, hess=ao_hess, deriv3=ao_d3, &
                               shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                               grad=ao_grad, hess=ao_hess, shell_mask=shell_mask, &
                               ao_offset=ao_offset, n_ao_out=n_sig)
         end if
         if (error%has_error()) return

         if (allocated(wg)) deallocate (wg)
         allocate (wg(nb))
         wg = ctx%grid%weights(g0:g1)

         ! (D chi) and (P chi): the two gemm partners everything below is
         ! assembled from, exactly `xc_hessian`'s `dmu` once per density.
         if (allocated(dmu)) deallocate (dmu, pmu)
         allocate (dmu(nb, n_sig), pmu(nb, n_sig))
         call pic_gemm(ao(1:nb, 1:n_sig), d_sig(1:n_sig, 1:n_sig), dmu, beta=0.0_dp)
         call pic_gemm(ao(1:nb, 1:n_sig), p_sig(1:n_sig, 1:n_sig), pmu, beta=0.0_dp)

         ! rho_P at every point, and for a GGA its gradient and the linearised
         ! sigma_P. Per block from the masked functions, against the *global*
         ! reference gradient -- the same mixed footing `xc_potential_gradient`
         ! stands on, so the two sides of the test truncate identically.
         if (allocated(rho_p)) deallocate (rho_p)
         allocate (rho_p(nb))
         do ig = 1, nb
            acc = 0.0_dp
            do i = 1, n_sig
               acc = acc + ao(ig, i)*pmu(ig, i)
            end do
            rho_p(ig) = acc
         end do

         if (gga) then
            if (allocated(dmu_d)) deallocate (dmu_d, pmu_d)
            allocate (dmu_d(nb, n_sig, 3), pmu_d(nb, n_sig, 3))
            do d = 1, 3
               call pic_gemm(ao_grad(1:nb, 1:n_sig, d), d_sig(1:n_sig, 1:n_sig), &
                             dmu_d(:, :, d), beta=0.0_dp)
               call pic_gemm(ao_grad(1:nb, 1:n_sig, d), p_sig(1:n_sig, 1:n_sig), &
                             pmu_d(:, :, d), beta=0.0_dp)
            end do
            if (allocated(p_grad)) deallocate (p_grad, sig_p)
            allocate (p_grad(nb, 3), sig_p(nb))
            do d = 1, 3
               do ig = 1, nb
                  acc = 0.0_dp
                  do i = 1, n_sig
                     acc = acc + ao_grad(ig, i, d)*pmu(ig, i)
                  end do
                  ! The two is the symmetry of P, exactly as in `eval_rho`.
                  p_grad(ig, d) = 2.0_dp*acc
               end do
            end do
            do ig = 1, nb
               gg = g0 + ig - 1
               sig_p(ig) = 2.0_dp*(rho_grad(gg, 1)*p_grad(ig, 1) &
                                   + rho_grad(gg, 2)*p_grad(ig, 2) &
                                   + rho_grad(gg, 3)*p_grad(ig, 3))
            end do
         end if

         ! The stand-in potentials: `f_xc` contracted with `(rho_P, sigma_P)`
         ! plays the role `v_xc` plays in `xc_hessian`'s potential terms.
         if (allocated(vrho_eff)) deallocate (vrho_eff)
         allocate (vrho_eff(nb))
         do ig = 1, nb
            vrho_eff(ig) = frr(g0 + ig - 1)*rho_p(ig)
         end do
         if (gga) then
            if (allocated(vsig_eff)) deallocate (vsig_eff)
            allocate (vsig_eff(nb))
            do ig = 1, nb
               gg = g0 + ig - 1
               vrho_eff(ig) = vrho_eff(ig) + frs(gg)*sig_p(ig)
               vsig_eff(ig) = frs(gg)*rho_p(ig) + fss(gg)*sig_p(ig)
            end do
         end if

         ! d rho/dA and d rho_P/dA: `xc_hessian`'s `dchi`, once per density.
         if (allocated(dchi)) deallocate (dchi, pchi)
         allocate (dchi(3*natm, nb), pchi(3*natm, nb))
         dchi = 0.0_dp
         pchi = 0.0_dp
         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               do ig = 1, nb
                  acc = 0.0_dp
                  accp = 0.0_dp
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     acc = acc + ao_grad(ig, i, a)*dmu(ig, i)
                     accp = accp + ao_grad(ig, i, a)*pmu(ig, i)
                  end do
                  dchi(3*(ia - 1) + a, ig) = -2.0_dp*acc
                  pchi(3*(ia - 1) + a, ig) = -2.0_dp*accp
               end do
            end do
         end do

         ! d(grad rho)/dA for both densities, then d sigma/dA and
         ! d sigma_P/dA. In the latter *both* factors of
         ! 2 grad rho . grad rho_P move -- forgetting the first product-rule
         ! half leaves an error confined to the sigma channel, which is what
         ! the B88 leg of the test exists to catch.
         if (gga) then
            if (allocated(dgrad)) deallocate (dgrad, pgrad, dsig, psig)
            allocate (dgrad(3*natm, 3, nb), pgrad(3*natm, 3, nb))
            allocate (dsig(3*natm, nb), psig(3*natm, nb))
            dgrad = 0.0_dp
            pgrad = 0.0_dp
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               do c = 1, 3
                  do d = 1, 3
                     do ig = 1, nb
                        acc = 0.0_dp
                        accp = 0.0_dp
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + ao_hess(ig, i, PAIR(c, d))*dmu(ig, i) &
                                 + ao_grad(ig, i, c)*dmu_d(ig, i, d)
                           accp = accp + ao_hess(ig, i, PAIR(c, d))*pmu(ig, i) &
                                  + ao_grad(ig, i, c)*pmu_d(ig, i, d)
                        end do
                        dgrad(3*(ia - 1) + c, d, ig) = -2.0_dp*acc
                        pgrad(3*(ia - 1) + c, d, ig) = -2.0_dp*accp
                     end do
                  end do
               end do
            end do
            do ig = 1, nb
               gg = g0 + ig - 1
               do p = 1, 3*natm
                  acc = 0.0_dp
                  accp = 0.0_dp
                  do d = 1, 3
                     acc = acc + rho_grad(gg, d)*dgrad(p, d, ig)
                     accp = accp + p_grad(ig, d)*dgrad(p, d, ig) &
                            + rho_grad(gg, d)*pgrad(p, d, ig)
                  end do
                  dsig(p, ig) = 2.0_dp*acc
                  psig(p, ig) = 2.0_dp*accp
               end do
            end do
         end if

         ! Groups one and two: the third derivatives against the reference's
         ! own first derivatives, and the kernel against the mixed ones. Same
         ! shape as `xc_hessian`'s kernel term with the ladder shifted a rung:
         ! where it reads `f`, this reads `g` contracted with `(rho_P,
         ! sigma_P)`, and the `f` that remains pairs a reference derivative
         ! with a `P` one.
         do ig = 1, nb
            gg = g0 + ig - 1
            wa = grrr(gg)*rho_p(ig)
            wb = 0.0_dp
            wc = 0.0_dp
            if (gga) then
               wa = wa + grrs(gg)*sig_p(ig)
               wb = grrs(gg)*rho_p(ig) + grss(gg)*sig_p(ig)
               wc = grss(gg)*rho_p(ig) + gsss(gg)*sig_p(ig)
            end if
            do p = 1, 3*natm
               do q = 1, 3*natm
                  acc = wa*dchi(p, ig)*dchi(q, ig) &
                        + frr(gg)*(dchi(p, ig)*pchi(q, ig) + pchi(p, ig)*dchi(q, ig))
                  if (gga) then
                     acc = acc &
                           + wb*(dchi(p, ig)*dsig(q, ig) + dsig(p, ig)*dchi(q, ig)) &
                           + wc*dsig(p, ig)*dsig(q, ig) &
                           + frs(gg)*(dchi(p, ig)*psig(q, ig) + psig(p, ig)*dchi(q, ig) &
                                      + dsig(p, ig)*pchi(q, ig) + pchi(p, ig)*dsig(q, ig)) &
                           + fss(gg)*(dsig(p, ig)*psig(q, ig) + psig(p, ig)*dsig(q, ig))
                  end if
                  hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) = &
                     hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) &
                     + wg(ig)*acc
               end do
            end do
         end do

         ! The `d(grad) . d(grad)` halves of the two d2-sigma terms: the
         ! kernel-weighted one over the reference's own gradient derivatives
         ! (group three), and the `v_sigma`-weighted mixed products from
         ! d2 sigma_P (group four).
         if (gga) then
            do ig = 1, nb
               gg = g0 + ig - 1
               do p = 1, 3*natm
                  do q = 1, 3*natm
                     acc = 0.0_dp
                     accp = 0.0_dp
                     do d = 1, 3
                        acc = acc + dgrad(p, d, ig)*dgrad(q, d, ig)
                        accp = accp + pgrad(p, d, ig)*dgrad(q, d, ig) &
                               + dgrad(p, d, ig)*pgrad(q, d, ig)
                     end do
                     hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) = &
                        hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) &
                        + 2.0_dp*wg(ig)*(vsig_eff(ig)*acc + vsigma(gg)*accp)
                  end do
               end do
            end do
         end if

         ! d2 rho against the kernel stand-in and d2 rho_P against the
         ! potential, off-diagonal piece: one derivative on each centre. Two
         ! passes over the same overlaps -- kernel weight against D, potential
         ! weight against P.
         if (allocated(zmat)) deallocate (zmat)
         allocate (zmat(n_sig, n_sig))
         do a = 1, 3
            do b = 1, 3
               call weighted_overlap(ao_grad(1:nb, 1:n_sig, a), ao_grad(1:nb, 1:n_sig, b), &
                                     wg*vrho_eff, zmat)
               call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                               natm, n_sig, a, b, 2.0_dp, diagonal=.false.)
               call weighted_overlap(ao_grad(1:nb, 1:n_sig, a), ao_grad(1:nb, 1:n_sig, b), &
                                     wg*vrho(g0:g1), zmat)
               call add_masked(hess, p_sig, zmat, c_offsets, c_counts, &
                               natm, n_sig, a, b, 2.0_dp, diagonal=.false.)
            end do
         end do

         ! The same two, diagonal piece: both derivatives on one function.
         do ih = 1, 6
            call hess_component(ih, a, b)
            call weighted_overlap(ao_hess(1:nb, 1:n_sig, ih), ao(1:nb, 1:n_sig), &
                                  wg*vrho_eff, zmat)
            call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                            natm, n_sig, a, b, 2.0_dp, diagonal=.true.)
            if (a /= b) call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                        natm, n_sig, b, a, 2.0_dp, diagonal=.true.)
            call weighted_overlap(ao_hess(1:nb, 1:n_sig, ih), ao(1:nb, 1:n_sig), &
                                  wg*vrho(g0:g1), zmat)
            call add_masked(hess, p_sig, zmat, c_offsets, c_counts, &
                            natm, n_sig, a, b, 2.0_dp, diagonal=.true.)
            if (a /= b) call add_masked(hess, p_sig, zmat, c_offsets, c_counts, &
                                        natm, n_sig, b, a, 2.0_dp, diagonal=.true.)
         end do

         ! The `grad . d2(grad)` halves, where the third basis-function
         ! derivatives live: three pairings where `xc_hessian` has one. The
         ! kernel stand-in with `grad rho` against D from d2 sigma; then, from
         ! d2 sigma_P, `v_sigma` with `grad rho_P` against D and `v_sigma`
         ! with `grad rho` against P -- the product rule on
         ! `grad rho . grad rho_P`, at second order this time.
         if (gga) then
            if (allocated(wsig)) deallocate (wsig)
            allocate (wsig(nb, 3))
            do d = 1, 3
               wsig(:, d) = 2.0_dp*wg*vsig_eff*rho_grad(g0:g1, d)
            end do
            call d2grad_sigma_term(hess, d_sig, wsig, ao, ao_grad, ao_hess, ao_d3, &
                                   c_offsets, c_counts, natm, n_sig, nb)
            do d = 1, 3
               wsig(:, d) = 2.0_dp*wg*vsigma(g0:g1)*p_grad(:, d)
            end do
            call d2grad_sigma_term(hess, d_sig, wsig, ao, ao_grad, ao_hess, ao_d3, &
                                   c_offsets, c_counts, natm, n_sig, nb)
            do d = 1, 3
               wsig(:, d) = 2.0_dp*wg*vsigma(g0:g1)*rho_grad(g0:g1, d)
            end do
            call d2grad_sigma_term(hess, p_sig, wsig, ao, ao_grad, ao_hess, ao_d3, &
                                   c_offsets, c_counts, natm, n_sig, nb)
         end if
      end do
   end subroutine xc_potential_hessian

   subroutine d2grad_sigma_term(hess, dens, wsig, ao, ao_grad, ao_hess, ao_d3, &
                                c_offsets, c_counts, natm, n_sig, nb)
      !! One `grad . d2(grad rho)` contraction: `sum_d wsig_d d2(grad rho)_d`
      !! over the density in `dens`, into the atom blocks of `hess`
      !!
      !! The four overlaps `xc_hessian` and `vv10_hessian` each spell out once:
      !! third derivative against the bare function plus second against first
      !! on the diagonal atom block -- both nuclear derivatives landing on one
      !! centre is what reaches a third position derivative -- and one
      !! derivative on each centre off it. `xc_potential_hessian` needs the
      !! same contraction three times with three different weight/density
      !! pairings, which is what made it a routine. Everything per point --
      !! quadrature weight, potential or kernel factor, whichever gradient is
      !! being dotted in -- travels in `wsig`.
      real(dp), intent(inout) :: hess(:, :, :, :)
      real(dp), intent(in) :: dens(:, :)
      real(dp), intent(in) :: wsig(:, :)   !! (nb, 3)
      real(dp), intent(in) :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), intent(in) :: ao_d3(:, :, :)
      integer, intent(in) :: c_offsets(:), c_counts(:)
      integer, intent(in) :: natm, n_sig, nb

      real(dp), allocatable :: zmat(:, :)
      integer :: c, d, e
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
      integer, parameter :: TRIP(3, 3, 3) = reshape( &
                            [1, 2, 3, 2, 4, 5, 3, 5, 6, &
                             2, 4, 5, 4, 7, 8, 5, 8, 9, &
                             3, 5, 6, 5, 8, 9, 6, 9, 10], [3, 3, 3])

      allocate (zmat(n_sig, n_sig))
      do d = 1, 3
         do c = 1, 3
            do e = 1, 3
               call weighted_overlap(ao_d3(1:nb, 1:n_sig, TRIP(c, e, d)), &
                                     ao(1:nb, 1:n_sig), wsig(:, d), zmat)
               call add_masked(hess, dens, zmat, c_offsets, c_counts, &
                               natm, n_sig, c, e, 2.0_dp, diagonal=.true.)
               call weighted_overlap(ao_hess(1:nb, 1:n_sig, PAIR(c, e)), &
                                     ao_grad(1:nb, 1:n_sig, d), wsig(:, d), zmat)
               call add_masked(hess, dens, zmat, c_offsets, c_counts, &
                               natm, n_sig, c, e, 2.0_dp, diagonal=.true.)
               call weighted_overlap(ao_hess(1:nb, 1:n_sig, PAIR(c, d)), &
                                     ao_grad(1:nb, 1:n_sig, e), wsig(:, d), zmat)
               call add_masked(hess, dens, zmat, c_offsets, c_counts, &
                               natm, n_sig, c, e, 2.0_dp, diagonal=.false.)
               call weighted_overlap(ao_grad(1:nb, 1:n_sig, c), &
                                     ao_hess(1:nb, 1:n_sig, PAIR(e, d)), wsig(:, d), zmat)
               call add_masked(hess, dens, zmat, c_offsets, c_counts, &
                               natm, n_sig, c, e, 2.0_dp, diagonal=.false.)
            end do
         end do
      end do
   end subroutine d2grad_sigma_term

   subroutine vv10_hessian(ctx, mol, density, hess, error)
      !! VV10's second derivative with the density matrix and both grids fixed
      !!
      !! The non-local counterpart of `xc_hessian`, on the NLC grid, and the
      !! second derivative of exactly what `vv10_gradient_fixed_grid` is the
      !! first derivative of. The same two-term shape holds, because `E_nl` is
      !! still a functional of rho and gamma alone:
      !!
      !!     d2 E / dA dB = int [ f_rho d2rho/dAdB + f_gamma d2gamma/dAdB ]
      !!                  + int int [ dA(rho,gamma) . f'' . dB(rho,gamma) ]
      !!
      !! The first term is the semilocal GGA contraction with `vrho` and
      !! `vsigma` swapped for VV10's -- Phase 1 pinned that its `f_rho` and
      !! `f_gamma` *are* our `vrho` and `vsigma` -- so it is `xc_hessian`'s
      !! potential machinery verbatim, minus the semilocal kernel. The second
      !! is where the non-locality lives: the kernel `f''` is a pair operator,
      !! so it is applied to all `3*natm` first derivatives at once by
      !! `vv10_hessian_kernel` and contracted here, following PySCF's
      !! `_get_enlc_deriv2` and the Liang, Feng, Liu & Head-Gordon paper it
      !! implements.
      !!
      !! The grid does not respond, deliberately, as in `xc_hessian` and as in
      !! PySCF's `grid_response=False` default -- which is why the check is a
      !! difference of the *fixed-grid* gradient and not the physical one.
      !! Restricted only, like the gradient upstream of it.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: hess(:, :, :, :)   !! (3, 3, natm, natm), in place
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), sigma(:), rho_grad(:, :)
      real(dp), allocatable :: rho_blk(:), rho_grad_blk(:, :)
      real(dp), allocatable :: exc(:), vrho(:), vsigma(:)
      real(dp), allocatable :: pu(:), pw(:), pa(:), pb(:), pc(:)
      real(dp), allocatable :: dodr(:), dodg(:), d2odr2(:), d2odg2(:), d2odrdg(:)
      real(dp), allocatable :: dkdr(:), d2kdr2(:)
      real(dp), allocatable :: drho_da(:, :), dgamma_da(:, :)
      real(dp), allocatable :: f_rho_t(:, :), f_gamma_t(:, :)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: ao_d3(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :)
      real(dp), allocatable :: dmu(:, :), dmu_d(:, :, :)
      real(dp), allocatable :: dgrad(:, :, :), wsig(:, :)
      real(dp), allocatable :: zmat(:, :), ymat(:, :), wg(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:), c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig
      integer :: ia, ja, a, b, i, j, ig, gg, p, q, c, d, e
      real(dp) :: acc
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
      integer, parameter :: TRIP(3, 3, 3) = reshape( &
                            [1, 2, 3, 2, 4, 5, 3, 5, 6, &
                             2, 4, 5, 4, 7, 8, 5, 8, 9, &
                             3, 5, 6, 5, 8, 9, 6, 9, 10], [3, 3, 3])

      if (error%has_error()) return
      if (ctx%nlc_b == 0.0_dp .and. ctx%nlc_c == 0.0_dp) return

      call ensure_nlc_grid(ctx, mol, error)
      if (error%has_error()) return
      npts = ctx%nlc_grid%n_points
      if (npts == 0) return

      natm = size(mol%coords, 2)
      nao = mol%nao

      ! Sweep one: rho, its gradient and sigma over the whole NLC grid, as the
      ! gradient does it -- the pair sums need every point before any output.
      allocate (rho(npts), sigma(npts), rho_grad(npts, 3))
      rho = 0.0_dp
      sigma = 0.0_dp
      rho_grad = 0.0_dp
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         if (allocated(rho_blk)) deallocate (rho_blk, rho_grad_blk)
         allocate (rho_blk(nb), rho_grad_blk(nb, 3))
         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, grad=ao_grad)
         if (error%has_error()) return
         call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=rho_grad_blk)
         do ig = 1, nb
            rho(g0 + ig - 1) = rho_blk(ig)
            do d = 1, 3
               rho_grad(g0 + ig - 1, d) = rho_grad_blk(ig, d)
            end do
            sigma(g0 + ig - 1) = rho_grad_blk(ig, 1)**2 + rho_grad_blk(ig, 2)**2 &
                                 + rho_grad_blk(ig, 3)**2
         end do
      end do

      ! One pair sweep for the potential and every kernel intermediate. The
      ! potential this returns *is* PySCF's `f_rho`/`f_gamma` -- Phase 1's test
      ! pins that identity -- so no rebuild happens here.
      allocate (exc(npts), vrho(npts), vsigma(npts))
      allocate (pu(npts), pw(npts), pa(npts), pb(npts), pc(npts))
      allocate (dodr(npts), dodg(npts), d2odr2(npts), d2odg2(npts), d2odrdg(npts))
      allocate (dkdr(npts), d2kdr2(npts))
      call vv10_nlc(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, rho, sigma, &
                    ctx%nlc_grid%coords, rho, sigma, ctx%nlc_grid%weights, &
                    exc, vrho, vsigma, &
                    hess_u=pu, hess_w=pw, hess_a=pa, hess_b=pb, hess_c=pc, &
                    domega_drho=dodr, domega_dgamma=dodg, &
                    d2omega_drho2=d2odr2, d2omega_dgamma2=d2odg2, &
                    d2omega_drho_dgamma=d2odrdg, &
                    dkappa_drho=dkdr, d2kappa_drho2=d2kdr2)

      ! The first derivatives of rho and gamma for all 3*natm perturbations,
      ! over the whole grid -- the kernel term consumes them globally, so they
      ! are stored rather than used and dropped per block as in `xc_hessian`.
      allocate (drho_da(3*natm, npts), dgamma_da(3*natm, npts))
      drho_da = 0.0_dp
      dgamma_da = 0.0_dp

      allocate (extents(mol%nbas))
      call shell_extents(mol, ctx%screen_tol, extents)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (c_offsets(natm), c_counts(natm))
      allocate (d_sig(nao, nao))

      ! Sweep two: term one -- the potential against the second derivatives of
      ! rho and gamma, `xc_hessian`'s machinery with its semilocal kernel
      ! removed -- and the storage of the first derivatives for term two.
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call block_significant_aos(mol, ctx%nlc_grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)
         if (n_sig == 0) cycle

         do jsig = 1, n_sig
            do isig = 1, n_sig
               d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
            end do
         end do

         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, &
                            grad=ao_grad, hess=ao_hess, deriv3=ao_d3, &
                            shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         if (error%has_error()) return

         if (allocated(wg)) deallocate (wg)
         allocate (wg(nb))
         wg = ctx%nlc_grid%weights(g0:g1)

         if (allocated(dmu)) deallocate (dmu)
         allocate (dmu(nb, n_sig))
         call pic_gemm(ao(1:nb, 1:n_sig), d_sig(1:n_sig, 1:n_sig), dmu, beta=0.0_dp)
         if (allocated(dmu_d)) deallocate (dmu_d)
         allocate (dmu_d(nb, n_sig, 3))
         do d = 1, 3
            call pic_gemm(ao_grad(1:nb, 1:n_sig, d), d_sig(1:n_sig, 1:n_sig), &
                          dmu_d(:, :, d), beta=0.0_dp)
         end do

         ! d rho / dA_c, straight into the global array.
         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               do ig = 1, nb
                  acc = 0.0_dp
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     acc = acc + ao_grad(ig, i, a)*dmu(ig, i)
                  end do
                  drho_da(3*(ia - 1) + a, g0 + ig - 1) = -2.0_dp*acc
               end do
            end do
         end do

         ! d(grad rho)/dA, and from it d gamma / dA.
         if (allocated(dgrad)) deallocate (dgrad)
         allocate (dgrad(3*natm, 3, nb))
         dgrad = 0.0_dp
         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do c = 1, 3
               do d = 1, 3
                  do ig = 1, nb
                     acc = 0.0_dp
                     do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                        acc = acc + ao_hess(ig, i, PAIR(c, d))*dmu(ig, i) &
                              + ao_grad(ig, i, c)*dmu_d(ig, i, d)
                     end do
                     dgrad(3*(ia - 1) + c, d, ig) = -2.0_dp*acc
                  end do
               end do
            end do
         end do
         do ig = 1, nb
            do p = 1, 3*natm
               acc = 0.0_dp
               do d = 1, 3
                  acc = acc + rho_grad(g0 + ig - 1, d)*dgrad(p, d, ig)
               end do
               dgamma_da(p, g0 + ig - 1) = 2.0_dp*acc
            end do
         end do

         ! Potential term against the second derivative of gamma, first half:
         ! the product of two first derivatives of grad rho.
         do ig = 1, nb
            do p = 1, 3*natm
               do q = 1, 3*natm
                  acc = 0.0_dp
                  do d = 1, 3
                     acc = acc + dgrad(p, d, ig)*dgrad(q, d, ig)
                  end do
                  hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) = &
                     hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) &
                     + 2.0_dp*wg(ig)*vsigma(g0 + ig - 1)*acc
               end do
            end do
         end do

         ! Potential term against d2rho, off-diagonal piece: one derivative on
         ! each centre.
         if (allocated(zmat)) deallocate (zmat)
         allocate (zmat(n_sig, n_sig))
         do a = 1, 3
            do b = 1, 3
               call weighted_overlap(ao_grad(1:nb, 1:n_sig, a), ao_grad(1:nb, 1:n_sig, b), &
                                     wg*vrho(g0:g1), zmat)
               do ja = 1, natm
                  if (c_counts(ja) == 0) cycle
                  do ia = 1, natm
                     if (c_counts(ia) == 0) cycle
                     acc = 0.0_dp
                     do j = c_offsets(ja) + 1, c_offsets(ja) + c_counts(ja)
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + d_sig(i, j)*zmat(i, j)
                        end do
                     end do
                     hess(a, b, ia, ja) = hess(a, b, ia, ja) + 2.0_dp*acc
                  end do
               end do
            end do
         end do

         ! Second half of the gamma term: `grad rho` against `d2(grad rho)`,
         ! third basis-function derivatives on the diagonal atom block.
         if (allocated(wsig)) deallocate (wsig)
         allocate (wsig(nb, 3))
         do d = 1, 3
            wsig(:, d) = 2.0_dp*wg*vsigma(g0:g1)*rho_grad(g0:g1, d)
         end do
         do d = 1, 3
            do c = 1, 3
               do e = 1, 3
                  call weighted_overlap(ao_d3(1:nb, 1:n_sig, TRIP(c, e, d)), &
                                        ao(1:nb, 1:n_sig), wsig(:, d), zmat)
                  call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                  natm, n_sig, c, e, 2.0_dp, diagonal=.true.)
                  call weighted_overlap(ao_hess(1:nb, 1:n_sig, PAIR(c, e)), &
                                        ao_grad(1:nb, 1:n_sig, d), wsig(:, d), zmat)
                  call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                  natm, n_sig, c, e, 2.0_dp, diagonal=.true.)
                  call weighted_overlap(ao_hess(1:nb, 1:n_sig, PAIR(c, d)), &
                                        ao_grad(1:nb, 1:n_sig, e), wsig(:, d), zmat)
                  call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                  natm, n_sig, c, e, 2.0_dp, diagonal=.false.)
                  call weighted_overlap(ao_grad(1:nb, 1:n_sig, c), &
                                        ao_hess(1:nb, 1:n_sig, PAIR(e, d)), wsig(:, d), zmat)
                  call add_masked(hess, d_sig, zmat, c_offsets, c_counts, &
                                  natm, n_sig, c, e, 2.0_dp, diagonal=.false.)
               end do
            end do
         end do

         ! Potential term against d2rho, diagonal piece: both derivatives on
         ! the same function.
         if (allocated(ymat)) deallocate (ymat)
         allocate (ymat(n_sig, n_sig))
         do p = 1, 6
            call hess_component(p, a, b)
            call weighted_overlap(ao_hess(1:nb, 1:n_sig, p), ao(1:nb, 1:n_sig), &
                                  wg*vrho(g0:g1), ymat)
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               acc = 0.0_dp
               do j = 1, n_sig
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     acc = acc + d_sig(i, j)*ymat(i, j)
                  end do
               end do
               hess(a, b, ia, ia) = hess(a, b, ia, ia) + 2.0_dp*acc
               if (a /= b) hess(b, a, ia, ia) = hess(b, a, ia, ia) + 2.0_dp*acc
            end do
         end do
      end do

      ! Term two: the kernel applied to every perturbation's first derivatives,
      ! then contracted against every other's -- one pair sweep for all 3*natm
      ! columns, then a weighted outer product over the grid. The inner
      ! quadrature weight lives inside the kernel; the outer one is applied
      ! here, exactly as PySCF splits them.
      allocate (f_rho_t(3*natm, npts), f_gamma_t(3*natm, npts))
      call vv10_hessian_kernel(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, &
                               rho, sigma, ctx%nlc_grid%weights, &
                               pu, pw, pa, pb, pc, dodr, dodg, dkdr, &
                               d2odr2, d2odg2, d2odrdg, d2kdr2, &
                               drho_da, dgamma_da, f_rho_t, f_gamma_t)
      do gg = 1, npts
         do q = 1, 3*natm
            do p = 1, 3*natm
               hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) = &
                  hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) &
                  + ctx%nlc_grid%weights(gg) &
                  *(drho_da(p, gg)*f_rho_t(q, gg) + dgamma_da(p, gg)*f_gamma_t(q, gg))
            end do
         end do
      end do
   end subroutine vv10_hessian

   subroutine xc_potential_deriv(ctx, mol, density, h1, error)
      !! The skeleton nuclear derivative of the exchange-correlation potential
      !!
      !! The response half of a Hessian contracts `dD/dR_B` against `dF/dR_A`,
      !! and for a Kohn-Sham reference `F` contains `V_xc`. Its derivative at
      !! fixed density matrix has two parts:
      !!
      !!     dV_uv/dA_c = int w v_rho d(chi_u chi_v)/dA_c
      !!                + int w f_rr chi_u chi_v (d rho/dA_c)
      !!
      !! the first from the basis functions moving and the second because the
      !! potential is itself a function of a density that moves with them.
      !!
      !! Leaving this out does not fail visibly. The coupled-perturbed solve
      !! converges, the Hessian comes back symmetric and plausible, and it is
      !! wrong by order one -- on water/STO-3G at LDA exchange the worst entry
      !! is out by 0.87 against a finite difference of the gradient, where the
      !! same comparison for Hartree-Fock agrees to 5e-6.
      !!
      !! A GGA adds the `sigma` channel to all of it: the potential itself gains
      !! `2 v_sigma grad rho . grad(chi_u chi_v)`, and every factor in that
      !! moves. Only second derivatives of the basis functions are needed here
      !! -- the third ones belong to the *energy's* second derivative, not the
      !! potential's first.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: h1(:, :, :, :)   !! (n_ao, n_ao, 3, natm), in place
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), rho_grad(:, :), vrho(:), vsigma(:)
      real(dp), allocatable :: frr(:), frs(:), fss(:)
      real(dp), allocatable :: tau(:), vtau(:), frt(:), fst(:), ftt(:)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), dmu(:, :), dchi(:, :)
      real(dp), allocatable :: wcol(:), blk(:, :), dmu_d(:, :, :)
      real(dp), allocatable :: dgrad(:, :, :), dsig(:, :)
      real(dp), allocatable :: gdot(:, :), dgdot(:, :), hgdot(:, :), dtau(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:), c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig, ia, a, i, j, ig, d
      real(dp) :: acc
      logical :: gga, mgga
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      if (error%has_error()) return

      natm = size(mol%coords, 2)
      nao = mol%nao
      npts = size(ctx%grid%weights)
      mgga = ctx%any_mgga
      gga = ctx%any_gga .or. mgga

      if (mgga) then
         call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                        frr, frs, fss, error, tau=tau, vtau=vtau, &
                                        frt=frt, fst=fst, ftt=ftt)
      else
         call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                        frr, frs, fss, error)
      end if
      if (error%has_error()) return

      allocate (extents(mol%nbas))
      call shell_extents(mol, ctx%screen_tol, extents)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (c_offsets(natm), c_counts(natm), d_sig(nao, nao))

      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)
         if (n_sig == 0) cycle
         do jsig = 1, n_sig
            do isig = 1, n_sig
               d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
            end do
         end do
         if (gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad, &
                               hess=ao_hess, shell_mask=shell_mask, ao_offset=ao_offset, &
                               n_ao_out=n_sig)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad, &
                               shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         end if
         if (error%has_error()) return

         if (allocated(dmu)) deallocate (dmu)
         allocate (dmu(nb, n_sig))
         call pic_gemm(ao(1:nb, 1:n_sig), d_sig(1:n_sig, 1:n_sig), dmu, beta=0.0_dp)

         ! d rho / dA_c on this block, as it is built in `xc_hessian`.
         if (allocated(dchi)) deallocate (dchi)
         allocate (dchi(3*natm, nb))
         dchi = 0.0_dp
         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               do ig = 1, nb
                  acc = 0.0_dp
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     acc = acc + ao_grad(ig, i, a)*dmu(ig, i)
                  end do
                  dchi(3*(ia - 1) + a, ig) = -2.0_dp*acc
               end do
            end do
         end do

         if (gga) then
            ! d(grad rho)/dA and d sigma/dA, as in `xc_hessian`, plus
            ! `gdot(g, u) = sum_d grad rho_d (d_d chi_u)`, which is the factor
            ! `grad rho . grad(chi_u chi_v)` reduces to once one index is fixed.
            if (allocated(dmu_d)) deallocate (dmu_d)
            allocate (dmu_d(nb, n_sig, 3))
            do d = 1, 3
               call pic_gemm(ao_grad(1:nb, 1:n_sig, d), d_sig(1:n_sig, 1:n_sig), &
                             dmu_d(:, :, d), beta=0.0_dp)
            end do
            if (allocated(dgrad)) deallocate (dgrad)
            if (allocated(dsig)) deallocate (dsig)
            if (allocated(gdot)) deallocate (gdot)
            allocate (dgrad(3*natm, 3, nb), dsig(3*natm, nb), gdot(nb, n_sig))
            dgrad = 0.0_dp
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               do a = 1, 3
                  do d = 1, 3
                     do ig = 1, nb
                        acc = 0.0_dp
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + ao_hess(ig, i, PAIR(a, d))*dmu(ig, i) &
                                 + ao_grad(ig, i, a)*dmu_d(ig, i, d)
                        end do
                        dgrad(3*(ia - 1) + a, d, ig) = -2.0_dp*acc
                     end do
                  end do
               end do
            end do
            do ig = 1, nb
               do i = 1, 3*natm
                  acc = 0.0_dp
                  do d = 1, 3
                     acc = acc + rho_grad(g0 + ig - 1, d)*dgrad(i, d, ig)
                  end do
                  dsig(i, ig) = 2.0_dp*acc
               end do
            end do
            gdot = 0.0_dp
            do d = 1, 3
               do i = 1, n_sig
                  do ig = 1, nb
                     gdot(ig, i) = gdot(ig, i) + rho_grad(g0 + ig - 1, d)*ao_grad(ig, i, d)
                  end do
               end do
            end do
         end if

         if (mgga) then
            if (allocated(dtau)) deallocate (dtau)
            allocate (dtau(3*natm, nb))
            dtau = 0.0_dp
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               do a = 1, 3
                  do ig = 1, nb
                     acc = 0.0_dp
                     do d = 1, 3
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + ao_hess(ig, i, PAIR(a, d))*dmu_d(ig, i, d)
                        end do
                     end do
                     dtau(3*(ia - 1) + a, ig) = -acc
                  end do
               end do
            end do
         end if

         if (allocated(blk)) deallocate (blk)
         if (allocated(wcol)) deallocate (wcol)
         allocate (blk(n_sig, n_sig), wcol(nb))
         if (gga) then
            if (allocated(dgdot)) deallocate (dgdot)
            if (allocated(hgdot)) deallocate (hgdot)
            allocate (dgdot(nb, n_sig), hgdot(nb, n_sig))
         end if

         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               ! Basis-function motion: minus the position derivative, and only
               ! on the rows this atom owns. The transpose supplies the ket half.
               wcol = ctx%grid%weights(g0:g1)*vrho(g0:g1)
               call weighted_overlap(ao_grad(1:nb, 1:n_sig, a), ao(1:nb, 1:n_sig), wcol, blk)
               do j = 1, n_sig
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     h1(ao_list(i), ao_list(j), a, ia) = &
                        h1(ao_list(i), ao_list(j), a, ia) - blk(i, j)
                     h1(ao_list(j), ao_list(i), a, ia) = &
                        h1(ao_list(j), ao_list(i), a, ia) - blk(i, j)
                  end do
               end do

               ! The potential following its own density. For a GGA `v_rho`
               ! depends on sigma as well, so the weight carries both channels.
               wcol = ctx%grid%weights(g0:g1)*frr(g0:g1)*dchi(3*(ia - 1) + a, 1:nb)
               if (gga) wcol = wcol + ctx%grid%weights(g0:g1)*frs(g0:g1) &
                               *dsig(3*(ia - 1) + a, 1:nb)
               if (mgga) wcol = wcol + ctx%grid%weights(g0:g1)*frt(g0:g1) &
                                *dtau(3*(ia - 1) + a, 1:nb)
               call weighted_overlap(ao(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
               call deposit_full(h1, ao_list, n_sig, a, ia, blk, 1.0_dp)

               if (gga) then
                  ! d v_sigma / dA, against `2 grad rho . grad(chi_u chi_v)`.
                  wcol = 2.0_dp*ctx%grid%weights(g0:g1) &
                         *(frs(g0:g1)*dchi(3*(ia - 1) + a, 1:nb) &
                           + fss(g0:g1)*dsig(3*(ia - 1) + a, 1:nb))
                  if (mgga) wcol = wcol + 2.0_dp*ctx%grid%weights(g0:g1)*fst(g0:g1) &
                                   *dtau(3*(ia - 1) + a, 1:nb)
                  call weighted_overlap(gdot(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
                  call deposit_full(h1, ao_list, n_sig, a, ia, blk, 1.0_dp)
                  call deposit_full(h1, ao_list, n_sig, a, ia, transpose(blk), 1.0_dp)

                  ! `grad rho` itself moving, against the same factor.
                  dgdot = 0.0_dp
                  do d = 1, 3
                     do i = 1, n_sig
                        do ig = 1, nb
                           dgdot(ig, i) = dgdot(ig, i) &
                                          + dgrad(3*(ia - 1) + a, d, ig)*ao_grad(ig, i, d)
                        end do
                     end do
                  end do
                  wcol = 2.0_dp*ctx%grid%weights(g0:g1)*vsigma(g0:g1)
                  call weighted_overlap(dgdot(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
                  call deposit_full(h1, ao_list, n_sig, a, ia, blk, 1.0_dp)
                  call deposit_full(h1, ao_list, n_sig, a, ia, transpose(blk), 1.0_dp)

                  ! The basis functions inside `grad(chi_u chi_v)` moving. Four
                  ! pieces, two on each index, and only the rows or columns this
                  ! atom owns.
                  hgdot = 0.0_dp
                  do d = 1, 3
                     do i = 1, n_sig
                        do ig = 1, nb
                           hgdot(ig, i) = hgdot(ig, i) &
                                          + rho_grad(g0 + ig - 1, d)*ao_hess(ig, i, PAIR(a, d))
                        end do
                     end do
                  end do
                  call weighted_overlap(hgdot(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
                  call deposit_rows(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                    c_offsets, c_counts)
                  call weighted_overlap(ao(1:nb, 1:n_sig), hgdot(1:nb, 1:n_sig), wcol, blk)
                  call deposit_cols(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                    c_offsets, c_counts)
                  call weighted_overlap(ao_grad(1:nb, 1:n_sig, a), gdot(1:nb, 1:n_sig), &
                                        wcol, blk)
                  call deposit_rows(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                    c_offsets, c_counts)
                  call weighted_overlap(gdot(1:nb, 1:n_sig), ao_grad(1:nb, 1:n_sig, a), &
                                        wcol, blk)
                  call deposit_cols(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                    c_offsets, c_counts)
               end if

               if (mgga) then
                  ! `v_tau` against half the gradient overlap of the two basis
                  ! functions. Its own derivative first, then the functions.
                  wcol = 0.5_dp*ctx%grid%weights(g0:g1) &
                         *(frt(g0:g1)*dchi(3*(ia - 1) + a, 1:nb) &
                           + fst(g0:g1)*dsig(3*(ia - 1) + a, 1:nb) &
                           + ftt(g0:g1)*dtau(3*(ia - 1) + a, 1:nb))
                  do d = 1, 3
                     call weighted_overlap(ao_grad(1:nb, 1:n_sig, d), &
                                           ao_grad(1:nb, 1:n_sig, d), wcol, blk)
                     call deposit_full(h1, ao_list, n_sig, a, ia, blk, 1.0_dp)
                  end do
                  wcol = 0.5_dp*ctx%grid%weights(g0:g1)*vtau(g0:g1)
                  do d = 1, 3
                     call weighted_overlap(ao_hess(1:nb, 1:n_sig, PAIR(a, d)), &
                                           ao_grad(1:nb, 1:n_sig, d), wcol, blk)
                     call deposit_rows(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                       c_offsets, c_counts)
                     call weighted_overlap(ao_grad(1:nb, 1:n_sig, d), &
                                           ao_hess(1:nb, 1:n_sig, PAIR(a, d)), wcol, blk)
                     call deposit_cols(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                       c_offsets, c_counts)
                  end do
               end if
            end do
         end do
      end do
   end subroutine xc_potential_deriv

   subroutine vv10_potential_deriv(ctx, mol, density, h1, error)
      !! The skeleton nuclear derivative of the VV10 potential matrix
      !!
      !! `xc_potential_deriv`'s non-local counterpart, on the NLC grid: VV10's
      !! contribution to `dF/dR` at fixed density matrix, for the response half
      !! of a Hessian to contract against `dD/dR`. The same two parts exist --
      !! the basis functions moving, and the potential being a functional of a
      !! density that moves with them -- but the second is non-local here. For
      !! a semilocal functional `d v_rho / dA` at a point is the kernel at that
      !! point against `d rho/dA` there; for VV10 a displaced density anywhere
      !! moves `f_rho` and `f_gamma` everywhere through the pair kernel, so the
      !! per-point weights `f_rr drho/dA + ...` become `vv10_hessian_kernel`'s
      !! `f_rho_t` and `f_gamma_t` -- one pair sweep for all `3*natm`
      !! perturbations, the same operator the energy's second derivative used.
      !! Which is also why those two deposits cannot skip atoms with no
      !! significant functions on a block: a distant atom's perturbation still
      !! changes the potential here, unlike every semilocal term.
      !!
      !! **Accumulates into `h1`**, like `xc_potential_deriv` and unlike
      !! `h1_contract` -- the caller zeroes.
      !!
      !! The grid does not respond, deliberately, matching `xc_potential_deriv`
      !! and PySCF's `grid_response=False` Hessian default; the orbital-response
      !! terms of PySCF's `_get_vnlc_deriv1` are what this transcribes, and the
      !! check is a difference of the fixed-grid potential for the reasons the
      !! Hessian tests set out. Restricted only, like everything upstream.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: h1(:, :, :, :)   !! (n_ao, n_ao, 3, natm), in place
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), sigma(:), rho_grad(:, :)
      real(dp), allocatable :: rho_blk(:), rho_grad_blk(:, :)
      real(dp), allocatable :: exc(:), vrho(:), vsigma(:)
      real(dp), allocatable :: pu(:), pw(:), pa(:), pb(:), pc(:)
      real(dp), allocatable :: dodr(:), dodg(:), d2odr2(:), d2odg2(:), d2odrdg(:)
      real(dp), allocatable :: dkdr(:), d2kdr2(:)
      real(dp), allocatable :: drho_da(:, :), dgamma_da(:, :)
      real(dp), allocatable :: f_rho_t(:, :), f_gamma_t(:, :)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :)
      real(dp), allocatable :: dmu(:, :), dmu_d(:, :, :), dgrad(:, :, :)
      real(dp), allocatable :: gdot(:, :), dgdot(:, :), hgdot(:, :)
      real(dp), allocatable :: wcol(:), blk(:, :), wg(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:), c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig
      integer :: ia, a, i, j, ig, d, p
      real(dp) :: acc
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      if (error%has_error()) return
      if (ctx%nlc_b == 0.0_dp .and. ctx%nlc_c == 0.0_dp) return

      call ensure_nlc_grid(ctx, mol, error)
      if (error%has_error()) return
      npts = ctx%nlc_grid%n_points
      if (npts == 0) return

      natm = size(mol%coords, 2)
      nao = mol%nao

      ! Sweep one: rho, its gradient and sigma over the whole NLC grid -- the
      ! pair sums need every point before any output, as in `vv10_hessian`.
      allocate (rho(npts), sigma(npts), rho_grad(npts, 3))
      rho = 0.0_dp
      sigma = 0.0_dp
      rho_grad = 0.0_dp
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         if (allocated(rho_blk)) deallocate (rho_blk, rho_grad_blk)
         allocate (rho_blk(nb), rho_grad_blk(nb, 3))
         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, grad=ao_grad)
         if (error%has_error()) return
         call eval_rho(ao, density, rho_blk, ao_grad=ao_grad, rho_grad=rho_grad_blk)
         do ig = 1, nb
            rho(g0 + ig - 1) = rho_blk(ig)
            do d = 1, 3
               rho_grad(g0 + ig - 1, d) = rho_grad_blk(ig, d)
            end do
            sigma(g0 + ig - 1) = rho_grad_blk(ig, 1)**2 + rho_grad_blk(ig, 2)**2 &
                                 + rho_grad_blk(ig, 3)**2
         end do
      end do

      ! One pair sweep for the potential -- which *is* PySCF's `f_rho` and
      ! `f_gamma`, the identity Phase 1's test pins -- and every kernel
      ! intermediate the perturbed potential needs.
      allocate (exc(npts), vrho(npts), vsigma(npts))
      allocate (pu(npts), pw(npts), pa(npts), pb(npts), pc(npts))
      allocate (dodr(npts), dodg(npts), d2odr2(npts), d2odg2(npts), d2odrdg(npts))
      allocate (dkdr(npts), d2kdr2(npts))
      call vv10_nlc(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, rho, sigma, &
                    ctx%nlc_grid%coords, rho, sigma, ctx%nlc_grid%weights, &
                    exc, vrho, vsigma, &
                    hess_u=pu, hess_w=pw, hess_a=pa, hess_b=pb, hess_c=pc, &
                    domega_drho=dodr, domega_dgamma=dodg, &
                    d2omega_drho2=d2odr2, d2omega_dgamma2=d2odg2, &
                    d2omega_drho_dgamma=d2odrdg, &
                    dkappa_drho=dkdr, d2kappa_drho2=d2kdr2)

      allocate (drho_da(3*natm, npts), dgamma_da(3*natm, npts))
      drho_da = 0.0_dp
      dgamma_da = 0.0_dp

      allocate (extents(mol%nbas))
      call shell_extents(mol, ctx%screen_tol, extents)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (c_offsets(natm), c_counts(natm))
      allocate (d_sig(nao, nao))

      ! Sweep two: the first derivatives of rho and gamma for the kernel, and
      ! every deposit that reads only the *unperturbed* potential -- the basis
      ! functions moving under `f_rho`, and the three ways the sigma-channel
      ! factor `2 f_gamma grad rho . grad(chi_u chi_v)` moves. These are
      ! `xc_potential_deriv`'s GGA terms verbatim with `v_rho`, `v_sigma`
      ! swapped for VV10's, and they are per-atom local, so the block's
      ! significant-atom mask applies to them.
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call block_significant_aos(mol, ctx%nlc_grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)
         if (n_sig == 0) cycle

         do jsig = 1, n_sig
            do isig = 1, n_sig
               d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
            end do
         end do

         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, &
                            grad=ao_grad, hess=ao_hess, &
                            shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         if (error%has_error()) return

         if (allocated(wg)) deallocate (wg)
         allocate (wg(nb))
         wg = ctx%nlc_grid%weights(g0:g1)

         if (allocated(dmu)) deallocate (dmu)
         allocate (dmu(nb, n_sig))
         call pic_gemm(ao(1:nb, 1:n_sig), d_sig(1:n_sig, 1:n_sig), dmu, beta=0.0_dp)
         if (allocated(dmu_d)) deallocate (dmu_d)
         allocate (dmu_d(nb, n_sig, 3))
         do d = 1, 3
            call pic_gemm(ao_grad(1:nb, 1:n_sig, d), d_sig(1:n_sig, 1:n_sig), &
                          dmu_d(:, :, d), beta=0.0_dp)
         end do

         ! d rho / dA_c, into the global array the kernel consumes.
         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               do ig = 1, nb
                  acc = 0.0_dp
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     acc = acc + ao_grad(ig, i, a)*dmu(ig, i)
                  end do
                  drho_da(3*(ia - 1) + a, g0 + ig - 1) = -2.0_dp*acc
               end do
            end do
         end do

         ! d(grad rho)/dA, and from it d gamma / dA.
         if (allocated(dgrad)) deallocate (dgrad)
         allocate (dgrad(3*natm, 3, nb))
         dgrad = 0.0_dp
         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               do d = 1, 3
                  do ig = 1, nb
                     acc = 0.0_dp
                     do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                        acc = acc + ao_hess(ig, i, PAIR(a, d))*dmu(ig, i) &
                              + ao_grad(ig, i, a)*dmu_d(ig, i, d)
                     end do
                     dgrad(3*(ia - 1) + a, d, ig) = -2.0_dp*acc
                  end do
               end do
            end do
         end do
         do ig = 1, nb
            do p = 1, 3*natm
               acc = 0.0_dp
               do d = 1, 3
                  acc = acc + rho_grad(g0 + ig - 1, d)*dgrad(p, d, ig)
               end do
               dgamma_da(p, g0 + ig - 1) = 2.0_dp*acc
            end do
         end do

         ! `gdot(g, u) = grad rho . grad chi_u`, shared by the sigma deposits.
         if (allocated(gdot)) deallocate (gdot)
         allocate (gdot(nb, n_sig))
         gdot = 0.0_dp
         do d = 1, 3
            do i = 1, n_sig
               do ig = 1, nb
                  gdot(ig, i) = gdot(ig, i) + rho_grad(g0 + ig - 1, d)*ao_grad(ig, i, d)
               end do
            end do
         end do

         if (allocated(blk)) deallocate (blk)
         if (allocated(wcol)) deallocate (wcol)
         if (allocated(dgdot)) deallocate (dgdot)
         if (allocated(hgdot)) deallocate (hgdot)
         allocate (blk(n_sig, n_sig), wcol(nb), dgdot(nb, n_sig), hgdot(nb, n_sig))

         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               ! Basis-function motion under `f_rho`: minus the position
               ! derivative, on the rows this atom owns, transpose for the ket.
               wcol = wg*vrho(g0:g1)
               call weighted_overlap(ao_grad(1:nb, 1:n_sig, a), ao(1:nb, 1:n_sig), wcol, blk)
               do j = 1, n_sig
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     h1(ao_list(i), ao_list(j), a, ia) = &
                        h1(ao_list(i), ao_list(j), a, ia) - blk(i, j)
                     h1(ao_list(j), ao_list(i), a, ia) = &
                        h1(ao_list(j), ao_list(i), a, ia) - blk(i, j)
                  end do
               end do

               ! `grad rho` itself moving, against `2 f_gamma grad(chi_u chi_v)`.
               dgdot = 0.0_dp
               do d = 1, 3
                  do i = 1, n_sig
                     do ig = 1, nb
                        dgdot(ig, i) = dgdot(ig, i) &
                                       + dgrad(3*(ia - 1) + a, d, ig)*ao_grad(ig, i, d)
                     end do
                  end do
               end do
               wcol = 2.0_dp*wg*vsigma(g0:g1)
               call weighted_overlap(dgdot(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
               call deposit_full(h1, ao_list, n_sig, a, ia, blk, 1.0_dp)
               call deposit_full(h1, ao_list, n_sig, a, ia, transpose(blk), 1.0_dp)

               ! The basis functions inside `grad(chi_u chi_v)` moving. Four
               ! pieces, two on each index, on the rows or columns this atom
               ! owns.
               hgdot = 0.0_dp
               do d = 1, 3
                  do i = 1, n_sig
                     do ig = 1, nb
                        hgdot(ig, i) = hgdot(ig, i) &
                                       + rho_grad(g0 + ig - 1, d)*ao_hess(ig, i, PAIR(a, d))
                     end do
                  end do
               end do
               call weighted_overlap(hgdot(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
               call deposit_rows(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                 c_offsets, c_counts)
               call weighted_overlap(ao(1:nb, 1:n_sig), hgdot(1:nb, 1:n_sig), wcol, blk)
               call deposit_cols(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                 c_offsets, c_counts)
               call weighted_overlap(ao_grad(1:nb, 1:n_sig, a), gdot(1:nb, 1:n_sig), &
                                     wcol, blk)
               call deposit_rows(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                 c_offsets, c_counts)
               call weighted_overlap(gdot(1:nb, 1:n_sig), ao_grad(1:nb, 1:n_sig, a), &
                                     wcol, blk)
               call deposit_cols(h1, ao_list, n_sig, a, ia, blk, -1.0_dp, &
                                 c_offsets, c_counts)
            end do
         end do
      end do

      ! The kernel: every perturbation's first derivatives in, the perturbed
      ! potentials `d f_rho / dA` and `d f_gamma / dA` out, in one pair sweep.
      ! The inner quadrature weight lives inside; the outer one is applied at
      ! the deposits below, exactly as PySCF weights `f_rho_t` only there.
      allocate (f_rho_t(3*natm, npts), f_gamma_t(3*natm, npts))
      call vv10_hessian_kernel(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, &
                               rho, sigma, ctx%nlc_grid%weights, &
                               pu, pw, pa, pb, pc, dodr, dodg, dkdr, &
                               d2odr2, d2odg2, d2odrdg, d2kdr2, &
                               drho_da, dgamma_da, f_rho_t, f_gamma_t)

      ! Sweep three: the perturbed potential against the *unmoved* basis pair,
      ! `w f_rho_t chi_u chi_v + 2 w f_gamma_t grad rho . grad(chi_u chi_v)`.
      ! **Every atom deposits on every block**: `f_rho_t` for atom A is nonzero
      ! wherever the grid is, whether or not A's functions reach this block --
      ! the pair kernel carried the perturbation there. The semilocal loop's
      ! `c_counts == 0` skip would silently drop exactly the long-range part.
      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1

         call block_significant_aos(mol, ctx%nlc_grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig)
         if (n_sig == 0) cycle

         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, error, &
                            grad=ao_grad, shell_mask=shell_mask, &
                            ao_offset=ao_offset, n_ao_out=n_sig)
         if (error%has_error()) return

         if (allocated(wg)) deallocate (wg)
         allocate (wg(nb))
         wg = ctx%nlc_grid%weights(g0:g1)

         if (allocated(gdot)) deallocate (gdot)
         allocate (gdot(nb, n_sig))
         gdot = 0.0_dp
         do d = 1, 3
            do i = 1, n_sig
               do ig = 1, nb
                  gdot(ig, i) = gdot(ig, i) + rho_grad(g0 + ig - 1, d)*ao_grad(ig, i, d)
               end do
            end do
         end do

         if (allocated(blk)) deallocate (blk)
         if (allocated(wcol)) deallocate (wcol)
         allocate (blk(n_sig, n_sig), wcol(nb))

         do ia = 1, natm
            do a = 1, 3
               p = 3*(ia - 1) + a
               wcol = wg*f_rho_t(p, g0:g1)
               call weighted_overlap(ao(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
               call deposit_full(h1, ao_list, n_sig, a, ia, blk, 1.0_dp)

               wcol = 2.0_dp*wg*f_gamma_t(p, g0:g1)
               call weighted_overlap(gdot(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
               call deposit_full(h1, ao_list, n_sig, a, ia, blk, 1.0_dp)
               call deposit_full(h1, ao_list, n_sig, a, ia, transpose(blk), 1.0_dp)
            end do
         end do
      end do
   end subroutine vv10_potential_deriv

   subroutine deposit_full(h1, ao_list, n_sig, a, ia, blk, scale)
      !! Add `scale * blk` into every element of one atom's perturbation block
      real(dp), intent(inout) :: h1(:, :, :, :)
      integer, intent(in) :: ao_list(:)
      integer, intent(in) :: n_sig, a, ia
      real(dp), intent(in) :: blk(:, :)
      real(dp), intent(in) :: scale
      integer :: i, j

      do j = 1, n_sig
         do i = 1, n_sig
            h1(ao_list(i), ao_list(j), a, ia) = h1(ao_list(i), ao_list(j), a, ia) &
                                                + scale*blk(i, j)
         end do
      end do
   end subroutine deposit_full

   subroutine deposit_rows(h1, ao_list, n_sig, a, ia, blk, scale, offsets, counts)
      !! As `deposit_full`, restricted to the rows atom `ia` owns
      real(dp), intent(inout) :: h1(:, :, :, :)
      integer, intent(in) :: ao_list(:)
      integer, intent(in) :: n_sig, a, ia
      real(dp), intent(in) :: blk(:, :)
      real(dp), intent(in) :: scale
      integer, intent(in) :: offsets(:), counts(:)

      integer :: i, j

      if (counts(ia) == 0) return
      do j = 1, n_sig
         do i = offsets(ia) + 1, offsets(ia) + counts(ia)
            h1(ao_list(i), ao_list(j), a, ia) = h1(ao_list(i), ao_list(j), a, ia) &
                                                + scale*blk(i, j)
         end do
      end do
   end subroutine deposit_rows

   subroutine deposit_cols(h1, ao_list, n_sig, a, ia, blk, scale, offsets, counts)
      !! As `deposit_full`, restricted to the columns atom `ia` owns
      real(dp), intent(inout) :: h1(:, :, :, :)
      integer, intent(in) :: ao_list(:)
      integer, intent(in) :: n_sig, a, ia
      real(dp), intent(in) :: blk(:, :)
      real(dp), intent(in) :: scale
      integer, intent(in) :: offsets(:), counts(:)

      integer :: i, j

      if (counts(ia) == 0) return
      do j = offsets(ia) + 1, offsets(ia) + counts(ia)
         do i = 1, n_sig
            h1(ao_list(i), ao_list(j), a, ia) = h1(ao_list(i), ao_list(j), a, ia) &
                                                + scale*blk(i, j)
         end do
      end do
   end subroutine deposit_cols

   subroutine xc_gradient_fixed_grid(ctx, mol, density, gradient, error)
      !! dE_xc/dR with the grid held fixed -- the first derivative `xc_hessian`
      !! is the second derivative of
      !!
      !! Not the physical exchange-correlation gradient: `xc_gradient` is, and it
      !! carries the grid-response terms this deliberately omits. This exists so
      !! the Hessian above can be differenced against its own first derivative,
      !! with the same approximation on both sides. Differencing the *physical*
      !! gradient instead would disagree by exactly the omitted term -- around
      !! 1e-4, small enough to look like a tolerance problem and large enough to
      !! hide a real error.
      type(xc_context_t), intent(inout) :: ctx
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(inout) :: gradient(:, :)   !! (3, natm), in place
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rho(:), rho_grad(:, :), vrho(:), vsigma(:)
      real(dp), allocatable :: frr(:), frs(:), fss(:)
      real(dp), allocatable :: tau(:), vtau(:), frt(:), fst(:), ftt(:)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), dmu(:, :), dmu_d(:, :, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:), c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig, ia, a, i, ig, d
      real(dp) :: acc
      logical :: gga, mgga
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      if (error%has_error()) return

      natm = size(mol%coords, 2)
      nao = mol%nao
      npts = size(ctx%grid%weights)

      mgga = ctx%any_mgga
      gga = ctx%any_gga .or. mgga
      if (mgga) then
         call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                        frr, frs, fss, error, tau=tau, vtau=vtau, &
                                        frt=frt, fst=fst, ftt=ftt)
      else
         call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                        frr, frs, fss, error)
      end if
      if (error%has_error()) return

      allocate (extents(mol%nbas))
      call shell_extents(mol, ctx%screen_tol, extents)
      allocate (shell_mask(mol%nbas), ao_offset(mol%nbas), ao_list(nao))
      allocate (c_offsets(natm), c_counts(natm), d_sig(nao, nao))

      do g0 = 1, npts, ctx%point_block
         g1 = min(g0 + ctx%point_block - 1, npts)
         nb = g1 - g0 + 1
         call block_significant_aos(mol, ctx%grid%coords(:, g0:g1), extents, &
                                    shell_mask, ao_list, ao_offset, n_sig, &
                                    atom_offsets=c_offsets, atom_counts=c_counts)
         if (n_sig == 0) cycle
         do jsig = 1, n_sig
            do isig = 1, n_sig
               d_sig(isig, jsig) = density(ao_list(isig), ao_list(jsig))
            end do
         end do
         if (gga) then
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad, &
                               hess=ao_hess, shell_mask=shell_mask, ao_offset=ao_offset, &
                               n_ao_out=n_sig)
         else
            call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad, &
                               shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         end if
         if (error%has_error()) return

         if (allocated(dmu)) deallocate (dmu)
         allocate (dmu(nb, n_sig))
         if (gga) then
            if (allocated(dmu_d)) deallocate (dmu_d)
            allocate (dmu_d(nb, n_sig, 3))
         end if
         call pic_gemm(ao(1:nb, 1:n_sig), d_sig(1:n_sig, 1:n_sig), dmu, beta=0.0_dp)

         do ia = 1, natm
            if (c_counts(ia) == 0) cycle
            do a = 1, 3
               acc = 0.0_dp
               do ig = 1, nb
                  do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                     acc = acc + ctx%grid%weights(g0 + ig - 1)*vrho(g0 + ig - 1) &
                           *ao_grad(ig, i, a)*dmu(ig, i)
                  end do
               end do
               gradient(a, ia) = gradient(a, ia) - 2.0_dp*acc
            end do
         end do

         ! The sigma channel: v_sigma against d sigma / dA, which is
         ! 2 grad rho . d(grad rho)/dA and costs one more derivative of every
         ! basis function than the rho term does.
         if (gga) then
            do d = 1, 3
               call pic_gemm(ao_grad(1:nb, 1:n_sig, d), d_sig(1:n_sig, 1:n_sig), &
                             dmu_d(:, :, d), beta=0.0_dp)
            end do
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               do a = 1, 3
                  acc = 0.0_dp
                  do d = 1, 3
                     do ig = 1, nb
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + 2.0_dp*ctx%grid%weights(g0 + ig - 1) &
                                 *vsigma(g0 + ig - 1)*rho_grad(g0 + ig - 1, d) &
                                 *(ao_hess(ig, i, PAIR(a, d))*dmu(ig, i) &
                                   + ao_grad(ig, i, a)*dmu_d(ig, i, d))
                        end do
                     end do
                  end do
                  gradient(a, ia) = gradient(a, ia) - 2.0_dp*acc
               end do
            end do
         end if

         if (mgga) then
            do ia = 1, natm
               if (c_counts(ia) == 0) cycle
               do a = 1, 3
                  acc = 0.0_dp
                  do d = 1, 3
                     do ig = 1, nb
                        do i = c_offsets(ia) + 1, c_offsets(ia) + c_counts(ia)
                           acc = acc + ctx%grid%weights(g0 + ig - 1)*vtau(g0 + ig - 1) &
                                 *ao_hess(ig, i, PAIR(a, d))*dmu_d(ig, i, d)
                        end do
                     end do
                  end do
                  gradient(a, ia) = gradient(a, ia) - acc
               end do
            end do
         end if
      end do
   end subroutine xc_gradient_fixed_grid

   subroutine add_masked(hess, dens, mat, offsets, counts, natm, n_sig, c, e, scale, diagonal)
      !! Add `scale * sum_{u in A, v in B} D_uv mat(u,v)` into `hess(c, e, A, B)`
      !!
      !! `diagonal` restricts the second index to every function rather than to
      !! atom B's, and puts the result on `hess(c, e, A, A)`. That is the case
      !! where both nuclear derivatives fell on the same centre, so there is no
      !! second atom to range over.
      real(dp), intent(inout) :: hess(:, :, :, :)
      real(dp), intent(in) :: dens(:, :), mat(:, :)
      integer, intent(in) :: offsets(:), counts(:)
      integer, intent(in) :: natm, n_sig, c, e
      real(dp), intent(in) :: scale
      logical, intent(in) :: diagonal

      integer :: ia, ja, i, j
      real(dp) :: acc

      do ia = 1, natm
         if (counts(ia) == 0) cycle
         if (diagonal) then
            acc = 0.0_dp
            do j = 1, n_sig
               do i = offsets(ia) + 1, offsets(ia) + counts(ia)
                  acc = acc + dens(i, j)*mat(i, j)
               end do
            end do
            hess(c, e, ia, ia) = hess(c, e, ia, ia) + scale*acc
         else
            do ja = 1, natm
               if (counts(ja) == 0) cycle
               acc = 0.0_dp
               do j = offsets(ja) + 1, offsets(ja) + counts(ja)
                  do i = offsets(ia) + 1, offsets(ia) + counts(ia)
                     acc = acc + dens(i, j)*mat(i, j)
                  end do
               end do
               hess(c, e, ia, ja) = hess(c, e, ia, ja) + scale*acc
            end do
         end if
      end do
   end subroutine add_masked

   subroutine weighted_overlap(left, right, w, out)
      !! `out(u,v) = sum_g w(g) left(g,u) right(g,v)`
      real(dp), intent(in) :: left(:, :), right(:, :), w(:)
      real(dp), intent(out) :: out(:, :)

      real(dp), allocatable :: scaled(:, :)
      integer :: ig, n

      n = size(left, 2)
      allocate (scaled(size(left, 1), n))
      do ig = 1, size(left, 1)
         scaled(ig, :) = w(ig)*left(ig, :)
      end do
      call pic_gemm(scaled, right, out, transa="T", beta=0.0_dp)
   end subroutine weighted_overlap

   pure subroutine hess_component(p, a, b)
      !! The packed second-derivative index `p` as its two Cartesian components
      !!
      !! `eval_ao_block` returns xx, xy, xz, yy, yz, zz -- the six unique
      !! entries of a symmetric tensor in `GTOval_sph_deriv2` order.
      integer, intent(in) :: p
      integer, intent(out) :: a, b

      select case (p)
      case (1)
         a = 1
         b = 1
      case (2)
         a = 1
         b = 2
      case (3)
         a = 1
         b = 3
      case (4)
         a = 2
         b = 2
      case (5)
         a = 2
         b = 3
      case default
         a = 3
         b = 3
      end select
   end subroutine hess_component

end module mqc_libcint_xc_hessian
