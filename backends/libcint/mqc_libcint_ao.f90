!! Basis function values on a grid, for DFT
module mqc_libcint_ao
   !! Evaluates every basis function of a molecule at a set of points.
   !!
   !! **libcint does not do this.** It computes integrals; there is no `GTOval` in
   !! it, and the function of that name belongs to PySCF's own `libcgto`. So this
   !! is the one piece of the DFT path that had to be written rather than called,
   !! and the whole difficulty is matching libcint's conventions exactly -- a basis
   !! function evaluated under a different-but-valid convention is a different
   !! basis, and the resulting SCF converges tidily onto a wrong number.
   !!
   !! Three conventions have to line up, and each is a way to be wrong quietly:
   !!
   !!   * **Normalisation is already in `env`.** `molecule_build` folds
   !!     `libcint_gto_norm` into every contraction coefficient, because libcint
   !!     expects that and does not apply it itself. Nothing here may normalise
   !!     again; doing so gives an overlap diagonal that is not one.
   !!   * **Cartesian component order** is libcint's: for l, the components run
   !!     `x^(l-i) y^(i-j) z^j` over `i = 0..l`, `j = 0..i`, which for d is
   !!     xx, xy, xz, yy, yz, zz.
   !!   * **The spherical transform must be libcint's**, coefficients and all. They
   !!     are in `mqc_libcint_ao_data`, transcribed from its own table, with the
   !!     s and p normalisation it keeps outside that table applied here.
   !!
   !! Which convention a molecule is in comes from `mol%cartesian`, the same flag
   !! `shell_dim` routes on, so the AO values cannot disagree with the AO count.
   !!
   !! Points are taken in blocks. A medium grid on a modest molecule is tens of
   !! thousands of points, and `n_points` by `n_ao` held whole is hundreds of
   !! megabytes at a real basis size -- the shape of mistake the coupled-cluster
   !! work already paid for once.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t, shell_dim
   use mqc_libcint_ao_data, only: C2S_LMAX, c2s_block, common_fac_sp
   use libcint_fortran, only: LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF
   implicit none
   private

   public :: eval_ao_block
   public :: shell_extents
   public :: block_significant_aos
   public :: eval_rho
   public :: max_ao_l

   !! Points per pass when a caller asks for a whole grid at once.
   integer, parameter, public :: AO_POINT_BLOCK = 512

   !! Unique components of a symmetric second-derivative tensor: xx, xy, xz,
   !! yy, yz, zz. Public because a caller indexing `hess` has to agree with the
   !! packing, and a bare 6 at both ends is how they stop agreeing.
   integer, parameter, public :: AO_HESS_COMP = 6
   integer, parameter, public :: AO_DERIV3_COMP = 10
      !! xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz -- the ten unique
      !! entries of a symmetric rank-three tensor, in `GTOval_sph_deriv3` order.

contains

   pure function max_ao_l(mol) result(l_max)
      !! Highest angular momentum in this molecule's basis
      type(libcint_molecule_t), intent(in) :: mol
      integer :: l_max

      integer :: ish

      l_max = 0
      do ish = 1, mol%nbas
         l_max = max(l_max, mol%bas(LIBCINT_ANG_OF, ish))
      end do
   end function max_ao_l

   subroutine shell_extents(mol, threshold, radius)
      !! Per shell, the distance past which every one of its functions is below
      !! `threshold`
      !!
      !! A contracted shell's value at distance r is a *sum* over primitives,
      !!
      !!     sum_p c_p r^l exp(-a_p r^2)
      !!
      !! so the envelope has to sum the coefficients, not take the largest of
      !! them. Every term is at most |c_p| r^l exp(-a_min r^2) since a_min is
      !! the exponent that survives longest, and summing those gives
      !!
      !!     (sum_p |c_p|) r^l exp(-a_min r^2)
      !!
      !! which is above the shell at every r. Taking `max_p |c_p|` instead --
      !! which this did first -- bounds one primitive rather than their sum and
      !! can be under the true value by up to a factor of `nprim`, worst on a
      !! heavily contracted s shell whose tightest primitives have not yet died
      !! at the cutoff. That made the stated threshold optimistic rather than
      !! the answer wrong, but "a bound" has to mean a bound.
      !! Solved by bisection rather than by the closed form: the closed form
      !! needs the Lambert W function for l > 0, and this is computed once per
      !! molecule and then reused for every block of every iteration, so a
      !! hundred bisection steps cost nothing worth avoiding.
      !!
      !! **A bound, not an estimate.** Every function of a shell is below the
      !! threshold beyond this radius, so dropping the shell outside it changes
      !! the quadrature by less than the threshold per point rather than by an
      !! amount nobody has characterised.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: threshold
      real(dp), allocatable, intent(out) :: radius(:)

      integer :: ish, l, nprim, nctr, ip, ic, exp_ptr, coef_ptr
      real(dp) :: a_min, c_max, c_sum, lo, hi, mid, val

      allocate (radius(mol%nbas))
      ! A non-positive threshold turns the screen off rather than being an
      ! error: it is how a run checks what the screening is worth, and how a
      ! disagreement gets bisected without rebuilding.
      if (threshold <= 0.0_dp) then
         radius = huge(1.0_dp)
         return
      end if
      do ish = 1, mol%nbas
         l = mol%bas(LIBCINT_ANG_OF, ish)
         nprim = mol%bas(LIBCINT_NPRIM_OF, ish)
         nctr = mol%bas(LIBCINT_NCTR_OF, ish)
         exp_ptr = mol%bas(LIBCINT_PTR_EXP, ish)
         coef_ptr = mol%bas(LIBCINT_PTR_COEFF, ish)

         a_min = huge(1.0_dp)
         do ip = 1, nprim
            a_min = min(a_min, mol%env(exp_ptr + ip))
         end do
         ! The largest contraction's total, since a general contraction's
         ! columns are separate functions and each has to be under the envelope.
         c_max = 0.0_dp
         do ic = 1, nctr
            c_sum = 0.0_dp
            do ip = 1, nprim
               c_sum = c_sum + abs(mol%env(coef_ptr + (ic - 1)*nprim + ip))
            end do
            c_max = max(c_max, c_sum)
         end do
         if (c_max <= 0.0_dp .or. a_min >= huge(1.0_dp)) then
            radius(ish) = 0.0_dp
            cycle
         end if

         ! Bracket first: r^l exp(-a r^2) is not monotone at small r for l > 0,
         ! so start past its maximum at r = sqrt(l/(2a)) and only then bisect.
         lo = sqrt(real(max(l, 1), dp)/(2.0_dp*a_min))
         hi = max(lo, 1.0_dp)
         do while (shell_bound(c_max, a_min, l, hi) > threshold .and. hi < 1.0e4_dp)
            hi = hi*2.0_dp
         end do
         do ip = 1, 80
            mid = 0.5_dp*(lo + hi)
            val = shell_bound(c_max, a_min, l, mid)
            if (val > threshold) then
               lo = mid
            else
               hi = mid
            end if
         end do
         radius(ish) = hi
      end do
   end subroutine shell_extents

   pure function shell_bound(c_max, a_min, l, r) result(v)
      !! |c| r^l exp(-a r^2), the envelope `shell_extents` inverts
      real(dp), intent(in) :: c_max, a_min, r
      integer, intent(in) :: l
      real(dp) :: v
      v = c_max*r**l*exp(-a_min*r*r)
   end function shell_bound

   subroutine block_significant_aos(mol, coords, radius, shell_mask, ao_list, &
                                    ao_offset, n_sig, atom_offsets, atom_counts)
      !! Which shells reach a block of points, and the AO indices they own
      !!
      !! The test is against the block's bounding sphere rather than each point:
      !! one distance per shell instead of `n_points`, and the bound stays a
      !! bound because a shell that fails it fails at every point in the block.
      !!
      !! Blocks that are spatially compact screen far harder than blocks that
      !! are not, which is the whole argument for keeping them small -- a block
      !! spanning the molecule reaches everything and screens nothing.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)     !! (3, n_points), Bohr
      real(dp), intent(in) :: radius(:)        !! from `shell_extents`
      logical, intent(out) :: shell_mask(:)    !! (nbas)
      integer, intent(out) :: ao_list(:)       !! (nao), first `n_sig` entries used
      integer, intent(out) :: ao_offset(:)
         !! (nbas) where each kept shell's functions start in the *compressed*
         !! numbering. Meaningless where `shell_mask` is false, and the reason
         !! `eval_ao_block` can write a compressed block without knowing which
         !! shells were dropped.
      integer, intent(out) :: n_sig
      integer, intent(out), optional :: atom_offsets(:), atom_counts(:)
         !! (natm) where each atom's kept functions begin in the compressed
         !! numbering, and how many it kept. Zero count for an atom that reaches
         !! this block with nothing.
         !!
         !! **These stay contiguous, and not by luck.** `molecule_build` packs
         !! shells with the loop over atoms outermost, so shells are strictly
         !! atom-ordered; the loop below walks them in that order handing out
         !! compressed indices as it goes, so an atom's kept functions are
         !! consecutive exactly as its full set was. `atom_ao_blocks` already
         !! relies on the same ordering for the unscreened case.
         !!
         !! The gradient is what needs them. `accumulate_channel` sums one
         !! atom's range onto that atom and every function onto the atom owning
         !! the point, with opposite signs, and those two cancel when the
         !! molecule translates. A range that straddled two atoms would put
         !! force on the wrong nucleus and leave a net force behind.

      integer :: ish, ig, npts, iatom, off_ao, ndim, i
      real(dp) :: centre(3)
      real(dp) :: span, d, dx, dy, dz

      npts = size(coords, 2)
      centre = 0.0_dp
      do ig = 1, npts
         centre = centre + coords(:, ig)
      end do
      centre = centre/real(max(npts, 1), dp)
      span = 0.0_dp
      do ig = 1, npts
         dx = coords(1, ig) - centre(1)
         dy = coords(2, ig) - centre(2)
         dz = coords(3, ig) - centre(3)
         span = max(span, sqrt(dx*dx + dy*dy + dz*dz))
      end do

      n_sig = 0
      off_ao = 0
      if (present(atom_offsets)) atom_offsets = 0
      if (present(atom_counts)) atom_counts = 0
      do ish = 1, mol%nbas
         ndim = shell_dim(mol%cartesian, ish - 1, mol%bas)
         iatom = mol%bas(LIBCINT_ATOM_OF, ish) + 1
         dx = mol%coords(1, iatom) - centre(1)
         dy = mol%coords(2, iatom) - centre(2)
         dz = mol%coords(3, iatom) - centre(3)
         d = sqrt(dx*dx + dy*dy + dz*dz)
         shell_mask(ish) = (d - span <= radius(ish))
         ao_offset(ish) = n_sig
         if (shell_mask(ish)) then
            if (present(atom_counts)) then
               ! The first kept shell of an atom fixes where its range starts.
               if (atom_counts(iatom) == 0 .and. present(atom_offsets)) then
                  atom_offsets(iatom) = n_sig
               end if
               atom_counts(iatom) = atom_counts(iatom) + ndim
            end if
            do i = 1, ndim
               n_sig = n_sig + 1
               ao_list(n_sig) = off_ao + i
            end do
         end if
         off_ao = off_ao + ndim
      end do
   end subroutine block_significant_aos

   subroutine eval_ao_block(mol, coords, ao, error, grad, hess, deriv3, shell_mask, &
                            ao_offset, n_ao_out)
      !! chi_mu(r) for every basis function at every supplied point
      !!
      !! `ao` comes back as (n_points, n_ao), which is the orientation the density
      !! assembly wants: rho at a point is a row contracted against the density
      !! matrix, so points vary slowest and a gemm over them reads contiguously.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)        !! (3, n_points), Bohr
      real(dp), allocatable, intent(out) :: ao(:, :)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: grad(:, :, :)
         !! d chi / d r, as (n_points, n_ao, 3). Asked for by GGA and above.
         !!
         !! Produced here rather than by a second routine because it shares
         !! everything expensive: the same exponentials, the same angular
         !! components, the same spherical transform. A separate `eval_ao_deriv1`
         !! would duplicate all three and have to be kept in step with them.
      logical, intent(in), optional :: shell_mask(:)
         !! (nbas) which shells reach these points, from `block_significant_aos`.
         !! Absent means all of them, which is the unscreened path every other
         !! caller still takes.
      integer, intent(in), optional :: ao_offset(:)
         !! (nbas) compressed offsets, required with `shell_mask`
      integer, intent(in), optional :: n_ao_out
         !! Leading dimension of the compressed block, required with `shell_mask`
      real(dp), allocatable, intent(out), optional :: deriv3(:, :, :)
         !! d3 chi / d r_i d r_j d r_k, as (n_points, n_ao, 10), packed
         !! xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz -- the order
         !! `GTOval_sph_deriv3` uses.
         !!
         !! A *GGA Hessian* needs these. The energy depends on grad rho, which
         !! already costs one position derivative of a basis function to form,
         !! and two nuclear derivatives on top of that land on the third. There
         !! is no cheaper identity; it is what the chain rule gives.
      real(dp), allocatable, intent(out), optional :: hess(:, :, :)
         !! d2 chi / d r_j d r_k, as (n_points, n_ao, 6), packed xx, xy, xz, yy,
         !! yz, zz -- the six unique components of a symmetric tensor, in the
         !! order PySCF's `GTOval_sph_deriv2` uses, so a comparison against it
         !! needs no reshuffling.
         !!
         !! Here for the same reason `grad` is: it shares the exponentials, the
         !! angular components and the transform with both of them. A GGA
         !! *gradient* needs these -- the energy depends on grad rho, and
         !! differentiating that with respect to a nuclear position differentiates
         !! the basis functions a second time.

      logical :: screened
      integer :: n_ao_here, shell_base
      integer :: n_points, ish, l, nprim, nctr, ic, ip, ig, off_ao
      integer :: n_cart, n_sph, n_here, comp, i, j, iatom
      integer :: exp_ptr, coef_ptr
      real(dp) :: dx, dy, dz, r2, radial, fac, dradial, d2radial
      real(dp) :: dr(3)
      real(dp) :: ex
      real(dp), allocatable :: cart(:), sph(:), trans(:, :)
      real(dp), allocatable :: dcart(:, :), dsph(:, :)
      real(dp), allocatable :: d2cart(:, :), d2sph(:, :)
      real(dp), allocatable :: d3cart(:, :), d3sph(:, :)
      real(dp) :: d3radial
      logical :: want_grad, want_hess, want_d3
      integer :: id, ih, px, py, pz

      ! The (j, k) each packed component stands for, and whether it is diagonal.
      ! The delta only enters the second derivative of the radial part -- see
      ! the assembly below -- so it is worth having rather than re-deriving.
      integer, parameter :: HESS_J(AO_HESS_COMP) = [1, 1, 1, 2, 2, 3]
      integer, parameter :: HESS_K(AO_HESS_COMP) = [1, 2, 3, 2, 3, 3]
      integer, parameter :: D3_I(AO_DERIV3_COMP) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
      integer, parameter :: D3_J(AO_DERIV3_COMP) = [1, 1, 1, 2, 2, 3, 2, 2, 3, 3]
      integer, parameter :: D3_K(AO_DERIV3_COMP) = [1, 2, 3, 2, 3, 3, 2, 3, 3, 3]
      integer, parameter :: HESS_OF(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])
         !! Which packed Hessian component a Cartesian pair (j,k) is.

      n_points = size(coords, 2)

      if (max_ao_l(mol) > C2S_LMAX) then
         call error%set(ERROR_VALIDATION, "AO evaluation: this basis has an angular "// &
                        "momentum above the transform table's range. Extend "// &
                        "mqc_libcint_ao_data rather than letting it through -- the "// &
                        "functions would otherwise be silently wrong.")
         return
      end if

      want_d3 = present(deriv3)
      ! Each order is built from the one below it, so asking for a third
      ! implies the second and the first.
      want_hess = present(hess) .or. want_d3
      ! The second derivative is built from the first, so asking for one
      ! implies the other. Computing grad silently and discarding it would be
      ! the alternative, and it is the same work.
      want_grad = present(grad) .or. want_hess
      screened = present(shell_mask)
      if (screened) then
         if (.not. (present(ao_offset) .and. present(n_ao_out))) then
            call error%set(ERROR_VALIDATION, "AO evaluation: a shell mask needs the "// &
                           "compressed offsets and width alongside it; passing one "// &
                           "without the others would silently write the wrong columns.")
            return
         end if
         n_ao_here = n_ao_out
      else
         n_ao_here = mol%nao
      end if

      allocate (ao(n_points, n_ao_here))
      ao = 0.0_dp
      if (present(grad)) then
         allocate (grad(n_points, n_ao_here, 3))
         grad = 0.0_dp
      end if
      if (present(hess)) then
         allocate (hess(n_points, n_ao_here, AO_HESS_COMP))
         hess = 0.0_dp
      end if
      if (want_d3) then
         allocate (deriv3(n_points, n_ao_here, AO_DERIV3_COMP))
         deriv3 = 0.0_dp
      end if

      !$omp parallel default(none) &
      !$omp    shared(mol, coords, ao, grad, hess, deriv3, n_points, want_grad, &
!$omp           want_hess, want_d3) &
      !$omp    shared(screened, shell_mask, ao_offset) &
      !$omp    private(ish, l, nprim, nctr, ic, ip, ig, off_ao, shell_base, &
      !$omp            n_cart, n_sph, &
      !$omp            n_here, comp, i, j, iatom, exp_ptr, coef_ptr, id, ih, &
      !$omp            px, py, pz, dr, ex, &
      !$omp            dx, dy, dz, r2, radial, dradial, d2radial, fac, cart, sph, trans, &
      !$omp            dcart, dsph, d2cart, d2sph, d3cart, d3sph, d3radial)
      allocate (cart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2), sph(2*C2S_LMAX + 1))
      allocate (trans(2*C2S_LMAX + 1, (C2S_LMAX + 1)*(C2S_LMAX + 2)/2))
      allocate (dcart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2, 3), dsph(2*C2S_LMAX + 1, 3))
      allocate (d2cart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2, AO_HESS_COMP), &
                d2sph(2*C2S_LMAX + 1, AO_HESS_COMP))
      allocate (d3cart((C2S_LMAX + 1)*(C2S_LMAX + 2)/2, AO_DERIV3_COMP), &
                d3sph(2*C2S_LMAX + 1, AO_DERIV3_COMP))

      do ish = 1, mol%nbas
         if (screened) then
            if (.not. shell_mask(ish)) cycle
            shell_base = ao_offset(ish)
         else
            shell_base = mol%shell_offset(ish)
         end if
         l = mol%bas(LIBCINT_ANG_OF, ish)
         nprim = mol%bas(LIBCINT_NPRIM_OF, ish)
         nctr = mol%bas(LIBCINT_NCTR_OF, ish)
         exp_ptr = mol%bas(LIBCINT_PTR_EXP, ish)
         coef_ptr = mol%bas(LIBCINT_PTR_COEFF, ish)
         iatom = mol%bas(LIBCINT_ATOM_OF, ish) + 1

         n_cart = (l + 1)*(l + 2)/2
         n_sph = 2*l + 1
         n_here = shell_dim(mol%cartesian, ish - 1, mol%bas)/nctr

         ! The factor libcint keeps outside its transform table, and it applies to
         ! *both* angular conventions -- `cint1e.c` multiplies every integral by
         ! `CINTcommon_fac_sp` of each shell's l without consulting cart-vs-sph.
         ! Gating it on the spherical path is wrong and shows up only in a
         ! Cartesian basis, where s and p come out too large by 1/0.2820948 and
         ! 1/0.4886025.
         fac = common_fac_sp(l)
         if (.not. mol%cartesian) call c2s_block(l, trans)

         !$omp do schedule(static)
         do ig = 1, n_points
            dx = coords(1, ig) - mol%coords(1, iatom)
            dy = coords(2, ig) - mol%coords(2, iatom)
            dz = coords(3, ig) - mol%coords(3, iatom)
            r2 = dx*dx + dy*dy + dz*dz
            dr = [dx, dy, dz]

            ! Angular part, in libcint's Cartesian order, and its gradient. The
            ! power rule term vanishes when the exponent is zero, which has to be
            ! a branch rather than an evaluation: pow(x, -1) is not what it means.
            comp = 0
            do i = 0, l
               do j = 0, i
                  comp = comp + 1
                  px = l - i
                  py = i - j
                  pz = j
                  cart(comp) = pow(dx, px)*pow(dy, py)*pow(dz, pz)
                  if (want_grad) then
                     dcart(comp, :) = 0.0_dp
                     if (px > 0) dcart(comp, 1) = real(px, dp) &
                                                  *pow(dx, px - 1)*pow(dy, py)*pow(dz, pz)
                     if (py > 0) dcart(comp, 2) = real(py, dp) &
                                                  *pow(dx, px)*pow(dy, py - 1)*pow(dz, pz)
                     if (pz > 0) dcart(comp, 3) = real(pz, dp) &
                                                  *pow(dx, px)*pow(dy, py)*pow(dz, pz - 1)
                  end if
                  if (want_hess) then
                     ! Same power rule applied twice. The guards are on the
                     ! *original* exponent rather than the reduced one, because
                     ! the prefactor already vanishes when it should: a first
                     ! power differentiated twice gives 1*0, and `pow` would
                     ! otherwise be asked for a negative exponent.
                     d2cart(comp, :) = 0.0_dp
                     if (px > 1) d2cart(comp, 1) = real(px*(px - 1), dp) &
                                                   *pow(dx, px - 2)*pow(dy, py)*pow(dz, pz)
                     if (px > 0 .and. py > 0) d2cart(comp, 2) = real(px*py, dp) &
                                                                *pow(dx, px - 1)*pow(dy, py - 1)*pow(dz, pz)
                     if (px > 0 .and. pz > 0) d2cart(comp, 3) = real(px*pz, dp) &
                                                                *pow(dx, px - 1)*pow(dy, py)*pow(dz, pz - 1)
                     if (py > 1) d2cart(comp, 4) = real(py*(py - 1), dp) &
                                                   *pow(dx, px)*pow(dy, py - 2)*pow(dz, pz)
                     if (py > 0 .and. pz > 0) d2cart(comp, 5) = real(py*pz, dp) &
                                                                *pow(dx, px)*pow(dy, py - 1)*pow(dz, pz - 1)
                     if (pz > 1) d2cart(comp, 6) = real(pz*(pz - 1), dp) &
                                                   *pow(dx, px)*pow(dy, py)*pow(dz, pz - 2)
                  end if
                  if (want_d3) then
                     ! The power rule a third time. `angular_third` keeps the
                     ! exponent bookkeeping in one place rather than writing ten
                     ! guarded expressions out, which is where a transposition
                     ! would hide.
                     do ih = 1, AO_DERIV3_COMP
                        d3cart(comp, ih) = angular_third(px, py, pz, dx, dy, dz, &
                                                         D3_I(ih), D3_J(ih), D3_K(ih))
                     end do
                  end if
               end do
            end do

            ! The transform is linear and does not depend on position, so it
            ! applies to the gradient exactly as it does to the value.
            if (.not. mol%cartesian) then
               do comp = 1, n_sph
                  sph(comp) = fac*sum(trans(comp, 1:n_cart)*cart(1:n_cart))
                  if (want_grad) then
                     do id = 1, 3
                        dsph(comp, id) = fac*sum(trans(comp, 1:n_cart)*dcart(1:n_cart, id))
                     end do
                  end if
                  if (want_hess) then
                     do ih = 1, AO_HESS_COMP
                        d2sph(comp, ih) = fac*sum(trans(comp, 1:n_cart)*d2cart(1:n_cart, ih))
                     end do
                  end if
                  if (want_d3) then
                     do ih = 1, AO_DERIV3_COMP
                        d3sph(comp, ih) = fac*sum(trans(comp, 1:n_cart)*d3cart(1:n_cart, ih))
                     end do
                  end if
               end do
            end if

            ! One contraction column at a time; libcint lays a shell's functions
            ! out with the contraction index outermost, and `molecule_build`
            ! preserved that, so the AO offset advances by n_here per column.
            do ic = 1, nctr
               radial = 0.0_dp
               dradial = 0.0_dp
               d2radial = 0.0_dp
               d3radial = 0.0_dp
               do ip = 1, nprim
                  ! d/dr_k exp(-a r^2) = -2 a r_k exp(-a r^2), so the -2a is
                  ! accumulated here and the r_k factor applied per component.
                  ! Differentiating again gives 4 a^2 r_j r_k - 2 a delta_jk;
                  ! the 4 a^2 is what accumulates below and the delta is applied
                  ! in the assembly, where it multiplies `dradial`.
                  ex = mol%env(coef_ptr + (ic - 1)*nprim + ip) &
                       *exp(-mol%env(exp_ptr + ip)*r2)
                  radial = radial + ex
                  if (want_grad) then
                     dradial = dradial - 2.0_dp*mol%env(exp_ptr + ip)*ex
                  end if
                  if (want_hess) then
                     d2radial = d2radial + 4.0_dp*mol%env(exp_ptr + ip) &
                                *mol%env(exp_ptr + ip)*ex
                  end if
                  if (want_d3) then
                     ! One more factor of -2a per derivative of the exponential.
                     d3radial = d3radial - 8.0_dp*mol%env(exp_ptr + ip)**3*ex
                  end if
               end do
               off_ao = shell_base + (ic - 1)*n_here
               if (mol%cartesian) then
                  do comp = 1, n_cart
                     ao(ig, off_ao + comp) = fac*radial*cart(comp)
                  end do
                  if (present(grad)) then
                     do comp = 1, n_cart
                        grad(ig, off_ao + comp, 1) = fac*(dcart(comp, 1)*radial &
                                                          + cart(comp)*dradial*dx)
                        grad(ig, off_ao + comp, 2) = fac*(dcart(comp, 2)*radial &
                                                          + cart(comp)*dradial*dy)
                        grad(ig, off_ao + comp, 3) = fac*(dcart(comp, 3)*radial &
                                                          + cart(comp)*dradial*dz)
                     end do
                  end if
                  if (want_d3) then
                     call assemble_third(d3cart, d2cart, dcart, cart, n_cart, &
                                         radial, dradial, d2radial, d3radial, dr, &
                                         fac, deriv3(ig, off_ao + 1:off_ao + n_cart, :))
                  end if
                  if (present(hess)) then
                     ! Product rule on chi = A(r) R(r^2), twice. The two cross
                     ! terms are not one term doubled: for an off-diagonal (j,k)
                     ! they differ, dA/dx_j pairing with dR/dx_k and vice versa.
                     do ih = 1, AO_HESS_COMP
                        i = HESS_J(ih)
                        j = HESS_K(ih)
                        hess(ig, off_ao + 1:off_ao + n_cart, ih) = &
                           fac*(d2cart(1:n_cart, ih)*radial &
                                + dcart(1:n_cart, i)*dradial*dr(j) &
                                + dcart(1:n_cart, j)*dradial*dr(i) &
                                + cart(1:n_cart)*d2radial*dr(i)*dr(j))
                        if (i == j) then
                           hess(ig, off_ao + 1:off_ao + n_cart, ih) = &
                              hess(ig, off_ao + 1:off_ao + n_cart, ih) &
                              + fac*cart(1:n_cart)*dradial
                        end if
                     end do
                  end if
               else
                  do comp = 1, n_sph
                     ao(ig, off_ao + comp) = radial*sph(comp)
                  end do
                  if (present(grad)) then
                     do comp = 1, n_sph
                        grad(ig, off_ao + comp, 1) = dsph(comp, 1)*radial &
                                                     + sph(comp)*dradial*dx
                        grad(ig, off_ao + comp, 2) = dsph(comp, 2)*radial &
                                                     + sph(comp)*dradial*dy
                        grad(ig, off_ao + comp, 3) = dsph(comp, 3)*radial &
                                                     + sph(comp)*dradial*dz
                     end do
                  end if
                  if (want_d3) then
                     call assemble_third(d3sph, d2sph, dsph, sph, n_sph, &
                                         radial, dradial, d2radial, d3radial, dr, &
                                         1.0_dp, deriv3(ig, off_ao + 1:off_ao + n_sph, :))
                  end if
                  if (present(hess)) then
                     ! As the Cartesian branch, with `fac` already inside sph,
                     ! dsph and d2sph rather than applied here.
                     do ih = 1, AO_HESS_COMP
                        i = HESS_J(ih)
                        j = HESS_K(ih)
                        hess(ig, off_ao + 1:off_ao + n_sph, ih) = &
                           d2sph(1:n_sph, ih)*radial &
                           + dsph(1:n_sph, i)*dradial*dr(j) &
                           + dsph(1:n_sph, j)*dradial*dr(i) &
                           + sph(1:n_sph)*d2radial*dr(i)*dr(j)
                        if (i == j) then
                           hess(ig, off_ao + 1:off_ao + n_sph, ih) = &
                              hess(ig, off_ao + 1:off_ao + n_sph, ih) &
                              + sph(1:n_sph)*dradial
                        end if
                     end do
                  end if
               end if
            end do
         end do
         !$omp end do
      end do

      deallocate (cart, sph, trans, dcart, dsph, d2cart, d2sph)
      !$omp end parallel
   end subroutine eval_ao_block

   pure subroutine assemble_third(a3, a2, a1, a0, n, r0, r1, r2, r3, dr, fac, out)
      !! The product rule on chi = A(r) R(r^2), three times
      !!
      !!   d3(AR) = A_ijk R
      !!          + A_ij R_k + A_ik R_j + A_jk R_i
      !!          + A_i R_jk + A_j R_ik + A_k R_ij
      !!          + A R_ijk
      !!
      !! with the radial factors following from R depending on r^2 alone:
      !!
      !!   R_i   = r1 x_i
      !!   R_ij  = r2 x_i x_j + r1 delta_ij
      !!   R_ijk = r3 x_i x_j x_k + r2 (delta_ij x_k + delta_ik x_j + delta_jk x_i)
      !!
      !! Shared by the Cartesian and spherical branches. The transform is linear
      !! and position-independent, so it commutes with every derivative -- which
      !! is why one routine serves both, with `fac` folded in on one side and
      !! already inside the angular arrays on the other.
      real(dp), intent(in) :: a3(:, :), a2(:, :), a1(:, :), a0(:)
      integer, intent(in) :: n
      real(dp), intent(in) :: r0, r1, r2, r3, fac
      real(dp), intent(in) :: dr(3)
      real(dp), intent(out) :: out(:, :)

      integer, parameter :: DI(AO_DERIV3_COMP) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
      integer, parameter :: DJ(AO_DERIV3_COMP) = [1, 1, 1, 2, 2, 3, 2, 2, 3, 3]
      integer, parameter :: DK(AO_DERIV3_COMP) = [1, 2, 3, 2, 3, 3, 2, 3, 3, 3]
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      integer :: ih, i, j, k
      real(dp) :: rk, rj, ri, rjk, rik, rij, rijk

      do ih = 1, AO_DERIV3_COMP
         i = DI(ih)
         j = DJ(ih)
         k = DK(ih)

         ri = r1*dr(i)
         rj = r1*dr(j)
         rk = r1*dr(k)

         rij = r2*dr(i)*dr(j)
         rik = r2*dr(i)*dr(k)
         rjk = r2*dr(j)*dr(k)
         if (i == j) rij = rij + r1
         if (i == k) rik = rik + r1
         if (j == k) rjk = rjk + r1

         rijk = r3*dr(i)*dr(j)*dr(k)
         if (i == j) rijk = rijk + r2*dr(k)
         if (i == k) rijk = rijk + r2*dr(j)
         if (j == k) rijk = rijk + r2*dr(i)

         out(1:n, ih) = fac*(a3(1:n, ih)*r0 &
                             + a2(1:n, PAIR(i, j))*rk &
                             + a2(1:n, PAIR(i, k))*rj &
                             + a2(1:n, PAIR(j, k))*ri &
                             + a1(1:n, i)*rjk &
                             + a1(1:n, j)*rik &
                             + a1(1:n, k)*rij &
                             + a0(1:n)*rijk)
      end do
   end subroutine assemble_third

   pure function angular_third(px, py, pz, dx, dy, dz, i, j, k) result(v)
      !! d3/dx_i dx_j dx_k of x^px y^py z^pz, at (dx, dy, dz)
      !!
      !! Written as exponent bookkeeping rather than ten guarded expressions:
      !! count how many derivatives fall on each axis, drop the term if any axis
      !! is differentiated more times than its exponent allows, and multiply the
      !! falling factorials. The guard is what keeps `pow` from being asked for
      !! a negative exponent, which is not what a vanishing prefactor means.
      integer, intent(in) :: px, py, pz, i, j, k
      real(dp), intent(in) :: dx, dy, dz
      real(dp) :: v

      integer :: n(3), p(3)
      integer :: m
      real(dp) :: coeff

      n = 0
      n(i) = n(i) + 1
      n(j) = n(j) + 1
      n(k) = n(k) + 1
      p = [px, py, pz]

      if (any(n > p)) then
         v = 0.0_dp
         return
      end if

      coeff = 1.0_dp
      do m = 1, 3
         coeff = coeff*falling(p(m), n(m))
      end do
      v = coeff*pow(dx, px - n(1))*pow(dy, py - n(2))*pow(dz, pz - n(3))
   end function angular_third

   pure function falling(p, n) result(f)
      !! p (p-1) ... (p-n+1), the factor n derivatives of x^p bring down
      integer, intent(in) :: p, n
      real(dp) :: f
      integer :: m

      f = 1.0_dp
      do m = 0, n - 1
         f = f*real(p - m, dp)
      end do
   end function falling

   subroutine eval_rho(ao, density, rho, ao_grad, rho_grad, tau)
      !! The electron density at every point, from the density matrix
      !!
      !!     rho(r) = sum_uv D_uv chi_u(r) chi_v(r)
      !!
      !! Through a gemm rather than a double loop over basis functions: form
      !! X = chi D once, then rho at a point is the row-wise dot of X with chi.
      !! The alternative is n_points n_ao^2 scalar work with the density matrix
      !! read in the worst possible order, which is the mistake the coupled-cluster
      !! ladder already paid for.
      !!
      !! `density` is whatever convention the caller's SCF uses -- the restricted
      !! path's D already carries two electrons per orbital, so this returns the
      !! total density for it and a spin density for an unrestricted D. Nothing
      !! here applies a factor, because doing so would be right for exactly one of
      !! those and silently wrong for the other.
      real(dp), intent(in) :: ao(:, :)        !! (n_points, n_ao)
      real(dp), intent(in) :: density(:, :)   !! (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: rho(:)
      real(dp), intent(in), optional :: ao_grad(:, :, :)
      real(dp), allocatable, intent(out), optional :: rho_grad(:, :)
         !! d rho / d r, as (n_points, 3). Both optionals or neither.
         !!
         !!     grad rho = 2 sum_uv D_uv chi_v grad chi_u
         !!
         !! The factor of two is the symmetry of D, not a spin factor, and it is
         !! the term that a GGA gets wrong by a few millihartree while converging
         !! perfectly.

      real(dp), allocatable, intent(out), optional :: tau(:)
         !! The kinetic energy density, for meta-GGA. Needs `ao_grad`.
         !!
         !!     tau(r) = 1/2 sum_d sum_uv D_uv d_d chi_u d_d chi_v
         !!
         !! The half is libxc's convention and PySCF's, checked against
         !! `eval_rho(..., xctype='MGGA')` rather than assumed -- the same quantity
         !! is defined without it elsewhere in the literature, and a factor of two
         !! in tau converges to a wrong energy like everything else here.

      real(dp), allocatable :: work(:, :), gwork(:, :)
      integer :: n_points, n_ao, ig, mu, id

      n_points = size(ao, 1)
      n_ao = size(ao, 2)

      allocate (work(n_points, n_ao), rho(n_points))
      call pic_gemm(ao, density, work)

      rho = 0.0_dp
      !$omp parallel do default(none) shared(rho, work, ao, n_points, n_ao) &
      !$omp    private(ig, mu) schedule(static)
      do ig = 1, n_points
         do mu = 1, n_ao
            rho(ig) = rho(ig) + work(ig, mu)*ao(ig, mu)
         end do
      end do
      !$omp end parallel do

      if (present(tau)) then
         allocate (tau(n_points))
         tau = 0.0_dp
         if (.not. present(ao_grad)) then
            tau = huge(1.0_dp)
         else
            allocate (gwork(n_points, n_ao))
            do id = 1, 3
               call pic_gemm(ao_grad(:, :, id), density, gwork)
               !$omp parallel do default(none) &
               !$omp    shared(tau, gwork, ao_grad, n_points, n_ao, id) &
               !$omp    private(ig, mu) schedule(static)
               do ig = 1, n_points
                  do mu = 1, n_ao
                     tau(ig) = tau(ig) + 0.5_dp*gwork(ig, mu)*ao_grad(ig, mu, id)
                  end do
               end do
               !$omp end parallel do
            end do
            deallocate (gwork)
         end if
      end if

      if (present(rho_grad)) then
         allocate (rho_grad(n_points, 3))
         rho_grad = 0.0_dp
         if (.not. present(ao_grad)) then
            ! Nothing to build it from; leaving zeros would be a silently
            ! LDA-shaped answer under a GGA's name.
            rho_grad = huge(1.0_dp)
         else
            !$omp parallel do default(none) &
            !$omp    shared(rho_grad, work, ao_grad, n_points, n_ao) &
            !$omp    private(ig, mu, id) schedule(static)
            do ig = 1, n_points
               do id = 1, 3
                  do mu = 1, n_ao
                     rho_grad(ig, id) = rho_grad(ig, id) &
                                        + 2.0_dp*work(ig, mu)*ao_grad(ig, mu, id)
                  end do
               end do
            end do
            !$omp end parallel do
         end if
      end if

      deallocate (work)
   end subroutine eval_rho

   pure function pow(x, n) result(v)
      !! x**n for small non-negative n, without the intrinsic's generality
      !!
      !! `x**0` must be one even when x is zero, which is exactly what happens at
      !! a grid point sitting on a nucleus -- and there are such points.
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: v

      integer :: k

      v = 1.0_dp
      do k = 1, n
         v = v*x
      end do
   end function pow

end module mqc_libcint_ao
