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
   use mqc_libcint_ao, only: eval_ao_block, shell_extents, block_significant_aos
   use mqc_libcint_xc, only: xc_context_t, xc_grid_kernel_quantities
   implicit none
   private

   public :: xc_hessian
   public :: xc_gradient_fixed_grid

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
      integer :: ia, ja, a, b, i, j, ig, p, q
      real(dp) :: acc

      if (error%has_error()) return

      if (ctx%any_gga .or. ctx%any_mgga) then
         call error%set(ERROR_VALIDATION, "An analytic Hessian is available for LDA "// &
                        "functionals only. A GGA needs third derivatives of the basis "// &
                        "functions, which the grid evaluator does not yet produce; use "// &
                        "the semi-numerical path for one.")
         return
      end if

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

         call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, &
                            grad=ao_grad, hess=ao_hess, shell_mask=shell_mask, &
                            ao_offset=ao_offset, n_ao_out=n_sig)
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

         ! Kernel term: f_rr contracted with the two first derivatives. An outer
         ! product over the block, which is one gemm rather than a loop nest.
         do ig = 1, nb
            acc = wg(ig)*frr(g0 + ig - 1)
            do p = 1, 3*natm
               do q = 1, 3*natm
                  hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) = &
                     hess(mod(p - 1, 3) + 1, mod(q - 1, 3) + 1, (p - 1)/3 + 1, (q - 1)/3 + 1) &
                     + acc*dchi(p, ig)*dchi(q, ig)
               end do
            end do
         end do

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
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: extents(:), d_sig(:, :), dmu(:, :)
      logical, allocatable :: shell_mask(:)
      integer, allocatable :: ao_offset(:), ao_list(:), c_offsets(:), c_counts(:)
      integer :: natm, nao, npts, g0, g1, nb, n_sig, isig, jsig, ia, a, i, ig
      real(dp) :: acc

      if (error%has_error()) return

      natm = size(mol%coords, 2)
      nao = mol%nao
      npts = size(ctx%grid%weights)

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
         call eval_ao_block(mol, ctx%grid%coords(:, g0:g1), ao, error, grad=ao_grad, &
                            shell_mask=shell_mask, ao_offset=ao_offset, n_ao_out=n_sig)
         if (error%has_error()) return

         if (allocated(dmu)) deallocate (dmu)
         allocate (dmu(nb, n_sig))
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
      end do
   end subroutine xc_gradient_fixed_grid

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
