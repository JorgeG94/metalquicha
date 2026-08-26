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
                             AO_DERIV3_COMP
   use mqc_libcint_xc, only: xc_context_t, xc_grid_kernel_quantities
   implicit none
   private

   public :: xc_hessian
   public :: xc_gradient_fixed_grid
   public :: xc_potential_deriv

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
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), dchi(:, :)
      real(dp), allocatable :: dmu(:, :), zmat(:, :), ymat(:, :), wg(:)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:)
      integer, allocatable :: c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig
      integer :: ia, ja, a, b, i, j, ig, p, q, c, d, e, ih
      real(dp) :: acc
      logical :: gga
      real(dp), allocatable :: ao_d3(:, :, :), dmu_d(:, :, :)
      real(dp), allocatable :: dgrad(:, :, :), dsig(:, :), wsig(:, :)
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
      integer, parameter :: TRIP(3, 3, 3) = reshape( &
                            [1, 2, 3, 2, 4, 5, 3, 5, 6, &
                             2, 4, 5, 4, 7, 8, 5, 8, 9, &
                             3, 5, 6, 5, 8, 9, 6, 9, 10], [3, 3, 3])
         !! Which packed third-derivative component a Cartesian triple is.

      if (error%has_error()) return

      if (ctx%any_mgga) then
         call error%set(ERROR_VALIDATION, "An analytic Hessian is available for LDA and "// &
                        "GGA functionals. A meta-GGA additionally depends on the kinetic "// &
                        "energy density, whose second derivative is not implemented; use "// &
                        "the semi-numerical path for one.")
         return
      end if
      gga = ctx%any_gga

      natm = size(mol%coords, 2)
      nao = mol%nao
      npts = size(ctx%grid%weights)

      ! rho, the potential and the kernel over the whole grid. The GGA channels
      ! come back too and are unused here, which the refusal above makes safe.
      call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                     frr, frs, fss, error)
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
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), dmu(:, :), dchi(:, :)
      real(dp), allocatable :: wcol(:), blk(:, :), dmu_d(:, :, :)
      real(dp), allocatable :: dgrad(:, :, :), dsig(:, :)
      real(dp), allocatable :: gdot(:, :), dgdot(:, :), hgdot(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:), c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig, ia, a, i, j, ig, d
      real(dp) :: acc
      logical :: gga
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      if (error%has_error()) return

      natm = size(mol%coords, 2)
      nao = mol%nao
      npts = size(ctx%grid%weights)
      gga = ctx%any_gga

      call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                     frr, frs, fss, error)
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
               call weighted_overlap(ao(1:nb, 1:n_sig), ao(1:nb, 1:n_sig), wcol, blk)
               call deposit_full(h1, ao_list, n_sig, a, ia, blk, 1.0_dp)

               if (gga) then
                  ! d v_sigma / dA, against `2 grad rho . grad(chi_u chi_v)`.
                  wcol = 2.0_dp*ctx%grid%weights(g0:g1) &
                         *(frs(g0:g1)*dchi(3*(ia - 1) + a, 1:nb) &
                           + fss(g0:g1)*dsig(3*(ia - 1) + a, 1:nb))
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
            end do
         end do
      end do
   end subroutine xc_potential_deriv

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
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :), ao_hess(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), dmu(:, :), dmu_d(:, :, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:), c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig, ia, a, i, ig, d
      real(dp) :: acc
      logical :: gga
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      if (error%has_error()) return

      natm = size(mol%coords, 2)
      nao = mol%nao
      npts = size(ctx%grid%weights)

      gga = ctx%any_gga
      call xc_grid_kernel_quantities(ctx, mol, density, rho, rho_grad, vrho, vsigma, &
                                     frr, frs, fss, error)
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
