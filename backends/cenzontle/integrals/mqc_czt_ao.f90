!! Basis function values on a grid, for DFT
module mqc_czt_ao
   !! Evaluates every basis function of a molecule at a set of points.
   !!
   !! **libcint does not do this.** It computes integrals; the `GTOval` of that
   !! name belongs to PySCF's own `libcgto`. So this is written rather than
   !! called, and the difficulty is matching libcint's conventions exactly -- a
   !! basis function evaluated under a different-but-valid convention is a
   !! different basis, and the SCF converges tidily onto a wrong number.
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
   !!     are in `mqc_czt_ao_data`, transcribed from its own table, with the
   !!     s and p normalisation it keeps outside that table applied here.
   !!
   !! Which convention a molecule is in comes from `mol%cartesian`, the same flag
   !! `shell_dim` routes on, so the AO values cannot disagree with the AO count.
   !!
   !! Points are taken in blocks: a medium grid on a modest molecule is tens of
   !! thousands of points, and `n_points` by `n_ao` held whole is hundreds of
   !! megabytes at a real basis size.
   use pic_types, only: dp
   use omp_lib, only: omp_get_num_threads, omp_get_thread_num
   use pic_blas_interfaces, only: pic_gemm
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, shell_dim
   use mqc_czt_ao_data, only: C2S_LMAX, c2s_block, common_fac_sp
   use libcint_fortran, only: LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, &
                              LIBCINT_NCTR_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF
   implicit none
   private

   public :: eval_ao_block
   public :: shell_extents
   public :: block_significant_aos
   public :: eval_rho
   public :: max_ao_l

   integer, parameter, public :: AO_POINT_BLOCK = 512
      !! Points per pass when a caller asks for a whole grid at once.

   integer, parameter, public :: AO_HESS_COMP = 6
      !! Unique components of a symmetric second-derivative tensor: xx, xy, xz,
      !! yy, yz, zz. Public because a caller indexing `hess` has to agree with
      !! the packing.
   integer, parameter, public :: AO_DERIV3_COMP = 10
      !! xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz -- the ten unique
      !! entries of a symmetric rank-three tensor, in `GTOval_sph_deriv3` order.

   integer, parameter :: AO_CHUNK = 128
      !! Points a shell is evaluated over at once: long enough to vectorise,
      !! short enough that the scratch stays in cache.
   real(dp), parameter :: EXP_CUTOFF = 50.0_dp
      !! A primitive is dropped at a point where a*r^2 exceeds this, i.e. where
      !! exp(-a r^2) < 2e-22. The callers' shell screening is 1e-12, so what is
      !! left out is ten orders of magnitude under what is already left out.
   integer, parameter :: N_DERIV_ALL = 1 + 3 + AO_HESS_COMP + AO_DERIV3_COMP
      !! Value, gradient, Hessian and third derivatives stacked: 20 components.
   integer, parameter :: DERIV_POWERS(3, N_DERIV_ALL) = reshape([ &
                                                                0, 0, 0, &
                                                                1, 0, 0, 0, 1, 0, 0, 0, 1, &
                                                                2, 0, 0, 1, 1, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 0, 2, &
                                                                3, 0, 0, 2, 1, 0, 2, 0, 1, 1, 2, 0, 1, 1, 1, 1, 0, 2, &
                                                                0, 3, 0, 0, 2, 1, 0, 1, 2, 0, 0, 3], [3, N_DERIV_ALL])
      !! How many times x, y and z are differentiated in each component of the
      !! stacked value, gradient, Hessian and third-derivative arrays, in the
      !! packing the public arrays use.

contains

   pure function max_ao_l(mol) result(l_max)
      !! Highest angular momentum in this molecule's basis
      type(czt_molecule_t), intent(in) :: mol
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
      !! which is above the shell at every r. Inverted by bisection, the closed
      !! form needing the Lambert W function for l > 0.
      !!
      !! **A bound, not an estimate.** Every function of a shell is below the
      !! threshold beyond this radius, so dropping the shell outside it changes
      !! the quadrature by less than the threshold per point.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: threshold
      real(dp), allocatable, intent(out) :: radius(:)

      integer :: ish, l, nprim, nctr, ip, ic, exp_ptr, coef_ptr
      real(dp) :: a_min, c_max, c_sum, lo, hi, mid, val

      allocate (radius(mol%nbas))
      ! A non-positive threshold turns the screen off rather than being an
      ! error, which is how a run checks what the screening is worth.
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
      !! A block spanning the molecule therefore reaches everything and screens
      !! nothing.
      type(czt_molecule_t), intent(in) :: mol
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
         !! **These stay contiguous.** `molecule_build` packs shells with the
         !! loop over atoms outermost, so shells are strictly atom-ordered and
         !! the loop below hands out compressed indices in that order. The
         !! gradient depends on it: a range straddling two atoms would put force
         !! on the wrong nucleus and leave a net force behind.

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
      !!
      !! Every element of every array is written, so none is zeroed first. Called
      !! from inside a parallel region it runs on the calling thread; called from
      !! serial code it splits the points over the team itself.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)        !! (3, n_points), Bohr
      real(dp), allocatable, intent(out) :: ao(:, :)
      type(error_t), intent(inout) :: error
      real(dp), allocatable, intent(out), optional :: grad(:, :, :)
         !! d chi / d r, as (n_points, n_ao, 3). Asked for by GGA and above.
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
         !! A *GGA Hessian* needs these: the energy depends on grad rho, which
         !! already costs one position derivative to form, and two nuclear
         !! derivatives on top of that land on the third.
      real(dp), allocatable, intent(out), optional :: hess(:, :, :)
         !! d2 chi / d r_j d r_k, as (n_points, n_ao, 6), packed xx, xy, xz, yy,
         !! yz, zz -- the six unique components of a symmetric tensor, in the
         !! order PySCF's `GTOval_sph_deriv2` uses, so a comparison against it
         !! needs no reshuffling. A GGA *gradient* needs these.

      logical :: screened, threaded
      integer :: n_ao_here, n_points, n_chunks, n_deriv, ish, kept

      n_points = size(coords, 2)

      if (max_ao_l(mol) > C2S_LMAX) then
         call error%set(ERROR_VALIDATION, "AO evaluation: this basis has an angular "// &
                        "momentum above the transform table's range. Extend "// &
                        "mqc_czt_ao_data rather than letting it through -- the "// &
                        "functions would otherwise be silently wrong.")
         return
      end if

      ! Each order is built from the one below it, so asking for a third
      ! implies the second and the second implies the first.
      n_deriv = 1
      if (present(grad)) n_deriv = 4
      if (present(hess)) n_deriv = 10
      if (present(deriv3)) n_deriv = 20
      screened = present(shell_mask)
      if (screened) then
         if (.not. (present(ao_offset) .and. present(n_ao_out))) then
            call error%set(ERROR_VALIDATION, "AO evaluation: a shell mask needs the "// &
                           "compressed offsets and width alongside it; passing one "// &
                           "without the others would silently write the wrong columns.")
            return
         end if
         ! Nothing below zeroes the outputs: every column is written by the
         ! shell that owns it. A width larger than the kept shells would leave
         ! columns holding whatever the allocator returned, so it is refused.
         kept = 0
         do ish = 1, mol%nbas
            if (shell_mask(ish)) kept = kept + shell_dim(mol%cartesian, ish - 1, mol%bas)
         end do
         if (n_ao_out /= kept) then
            call error%set(ERROR_VALIDATION, "AO evaluation: n_ao_out does not equal the "// &
                           "number of functions in the masked shells, so some output "// &
                           "columns would never be written.")
            return
         end if
         n_ao_here = n_ao_out
      else
         n_ao_here = mol%nao
      end if

      allocate (ao(n_points, n_ao_here))
      if (present(grad)) allocate (grad(n_points, n_ao_here, 3))
      if (present(hess)) allocate (hess(n_points, n_ao_here, AO_HESS_COMP))
      if (present(deriv3)) allocate (deriv3(n_points, n_ao_here, AO_DERIV3_COMP))

      ! The callers sit inside their own parallel loop over blocks of points, so
      ! the common case is a team already running and nothing to open here. A
      ! caller from serial code gets the team, over chunks of points.
      n_chunks = (n_points + AO_CHUNK - 1)/AO_CHUNK
      threaded = n_chunks > 1 .and. omp_get_num_threads() == 1
      !$omp parallel default(shared) if(threaded)
      call eval_chunks(mol, coords, n_chunks, omp_get_thread_num(), omp_get_num_threads(), &
                       n_deriv, screened, shell_mask, ao_offset, ao, grad, hess, deriv3)
      !$omp end parallel
   end subroutine eval_ao_block

   subroutine eval_chunks(mol, coords, n_chunks, first, stride, n_deriv, screened, &
                          shell_mask, ao_offset, ao, grad, hess, deriv3)
      !! Chunks `first+1, first+1+stride, ...` of an `eval_ao_block` call
      !!
      !! Holds the scratch for those chunks, sized once from the molecule's
      !! highest l and reused for every chunk and shell.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: coords(:, :)
      integer, intent(in) :: n_chunks, first, stride
      integer, intent(in) :: n_deriv
         !! 1, 4, 10 or 20: how many of value, gradient, Hessian, third
         !! derivatives are built, each including the ones below it
      logical, intent(in) :: screened
      logical, intent(in), optional :: shell_mask(:)
      integer, intent(in), optional :: ao_offset(:)
      real(dp), intent(inout) :: ao(:, :)
      real(dp), intent(inout), optional :: grad(:, :, :)
      real(dp), intent(inout), optional :: hess(:, :, :)
      real(dp), intent(inout), optional :: deriv3(:, :, :)

      integer :: l_max, n_cart_max, n_sph_max, n_rad
      integer :: ichunk, p0, np, ish, l, nprim, nctr, ic, ip, iatom, k, p
      integer :: n_cart, n_sph, n_fun, n_here, shell_base, off_ao, comp, s, d
      integer :: exp_ptr, coef_ptr, i, j, n_sel, px, py, pz, e
      integer :: n_pow(3)
      real(dp) :: ax, ay, az, expo, coef, ar2, t, p_coef
      real(dp), allocatable :: dd(:, :), r2(:), pw(:, :, :), rad(:, :), arg(:), ev(:)
      real(dp), allocatable :: cart(:, :, :), sph(:, :, :), trans(:, :, :), tmp(:, :)
      integer, allocatable :: idx(:)

      l_max = max_ao_l(mol)
      n_cart_max = (l_max + 1)*(l_max + 2)/2
      n_sph_max = 2*l_max + 1
      n_rad = 1
      if (n_deriv >= 4) n_rad = 2
      if (n_deriv >= 10) n_rad = 3
      if (n_deriv >= 20) n_rad = 4

      allocate (dd(AO_CHUNK, 3), r2(AO_CHUNK), pw(AO_CHUNK, 0:l_max, 3), rad(AO_CHUNK, 4))
      allocate (arg(AO_CHUNK), ev(AO_CHUNK), idx(AO_CHUNK))
      allocate (cart(AO_CHUNK, n_cart_max, n_deriv))
      if (.not. mol%cartesian) then
         allocate (sph(AO_CHUNK, n_sph_max, n_deriv))
         allocate (trans(n_sph_max, n_cart_max, 0:l_max))
         trans = 0.0_dp
         do l = 0, l_max
            allocate (tmp(2*l + 1, (l + 1)*(l + 2)/2))
            call c2s_block(l, tmp)
            trans(1:2*l + 1, 1:(l + 1)*(l + 2)/2, l) = tmp
            deallocate (tmp)
         end do
      end if

      do ichunk = first + 1, n_chunks, stride
         p0 = (ichunk - 1)*AO_CHUNK
         np = min(AO_CHUNK, size(coords, 2) - p0)

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
            n_fun = n_cart
            if (.not. mol%cartesian) n_fun = n_sph

            ax = mol%coords(1, iatom)
            ay = mol%coords(2, iatom)
            az = mol%coords(3, iatom)
            do p = 1, np
               dd(p, 1) = coords(1, p0 + p) - ax
               dd(p, 2) = coords(2, p0 + p) - ay
               dd(p, 3) = coords(3, p0 + p) - az
               r2(p) = dd(p, 1)*dd(p, 1) + dd(p, 2)*dd(p, 2) + dd(p, 3)*dd(p, 3)
            end do

            ! x^0 .. x^l by repeated multiplication. The zeroth power is set,
            ! not computed, so a point sitting on the nucleus still gets one.
            do d = 1, 3
               pw(1:np, 0, d) = 1.0_dp
               do e = 1, l
                  pw(1:np, e, d) = pw(1:np, e - 1, d)*dd(1:np, d)
               end do
            end do

            ! Angular part in libcint's Cartesian order, with every derivative
            ! taken from the same power tables. `common_fac_sp` is folded in
            ! here; it applies to both angular conventions, and the transform
            ! below is linear.
            comp = 0
            do i = 0, l
               do j = 0, i
                  comp = comp + 1
                  px = l - i
                  py = i - j
                  pz = j
                  do d = 1, n_deriv
                     n_pow = DERIV_POWERS(:, d)
                     if (n_pow(1) > px .or. n_pow(2) > py .or. n_pow(3) > pz) then
                        cart(1:np, comp, d) = 0.0_dp
                        cycle
                     end if
                     p_coef = common_fac_sp(l)*falling(px, n_pow(1))*falling(py, n_pow(2)) &
                              *falling(pz, n_pow(3))
                     cart(1:np, comp, d) = p_coef*pw(1:np, px - n_pow(1), 1) &
                                           *pw(1:np, py - n_pow(2), 2)*pw(1:np, pz - n_pow(3), 3)
                  end do
               end do
            end do

            ! The transform is linear and does not depend on position, so it
            ! applies to every derivative as it does to the value. Most of its
            ! coefficients are zero and are skipped.
            if (.not. mol%cartesian) then
               do d = 1, n_deriv
                  do s = 1, n_sph
                     sph(1:np, s, d) = 0.0_dp
                     do comp = 1, n_cart
                        t = trans(s, comp, l)
                        if (t == 0.0_dp) cycle
                        sph(1:np, s, d) = sph(1:np, s, d) + t*cart(1:np, comp, d)
                     end do
                  end do
               end do
            end if

            ! One contraction column at a time: libcint lays a shell's
            ! functions out with the contraction index outermost, so the AO
            ! offset advances by `n_here` per column.
            do ic = 1, nctr
               rad(1:np, 1:n_rad) = 0.0_dp
               do ip = 1, nprim
                  expo = mol%env(exp_ptr + ip)
                  coef = mol%env(coef_ptr + (ic - 1)*nprim + ip)
                  ! Only where the exponential is above e^-EXP_CUTOFF; gathered
                  ! so that the `exp` runs over a dense list.
                  n_sel = 0
                  do p = 1, np
                     ar2 = expo*r2(p)
                     if (ar2 <= EXP_CUTOFF) then
                        n_sel = n_sel + 1
                        idx(n_sel) = p
                        arg(n_sel) = -ar2
                     end if
                  end do
                  !$omp simd
                  do k = 1, n_sel
                     ev(k) = coef*exp(arg(k))
                  end do
                  ! d/dr_k exp(-a r^2) = -2 a r_k exp(-a r^2), so the -2a is
                  ! accumulated here and the r_k factor applied per component.
                  ! Differentiating again gives 4 a^2 r_j r_k - 2 a delta_jk;
                  ! the 4 a^2 is what accumulates below and the delta is applied
                  ! in the assembly, where it multiplies `rad(:, 2)`.
                  do k = 1, n_sel
                     rad(idx(k), 1) = rad(idx(k), 1) + ev(k)
                  end do
                  if (n_rad >= 2) then
                     do k = 1, n_sel
                        rad(idx(k), 2) = rad(idx(k), 2) - 2.0_dp*expo*ev(k)
                     end do
                  end if
                  if (n_rad >= 3) then
                     do k = 1, n_sel
                        rad(idx(k), 3) = rad(idx(k), 3) + 4.0_dp*expo*expo*ev(k)
                     end do
                  end if
                  if (n_rad >= 4) then
                     ! One more factor of -2a per derivative of the exponential.
                     do k = 1, n_sel
                        rad(idx(k), 4) = rad(idx(k), 4) - 8.0_dp*expo**3*ev(k)
                     end do
                  end if
               end do

               off_ao = shell_base + (ic - 1)*n_here
               if (mol%cartesian) then
                  call store_shell(cart, n_fun, np, p0, off_ao, rad, dd, ao, grad, hess, deriv3)
               else
                  call store_shell(sph, n_fun, np, p0, off_ao, rad, dd, ao, grad, hess, deriv3)
               end if
            end do
         end do
      end do
   end subroutine eval_chunks

   subroutine store_shell(a, n_fun, np, p0, off_ao, rad, dd, ao, grad, hess, deriv3)
      !! Angular times radial, by the product rule, into the output columns
      !!
      !! chi = A(r) R(r^2). The radial factors follow from R depending on r^2
      !! alone:
      !!
      !!   R_i   = r1 x_i
      !!   R_ij  = r2 x_i x_j + r1 delta_ij
      !!   R_ijk = r3 x_i x_j x_k + r2 (delta_ij x_k + delta_ik x_j + delta_jk x_i)
      !!
      !! and the angular derivatives are in `a`, packed as `DERIV_POWERS`. The
      !! two cross terms of a second derivative are not one term doubled: for an
      !! off-diagonal (j,k) they differ, dA/dx_j pairing with dR/dx_k and vice
      !! versa.
      real(dp), intent(in) :: a(:, :, :)
         !! (chunk, n_fun, n_deriv) angular parts and their derivatives, with the
         !! shell's normalisation inside
      integer, intent(in) :: n_fun, np, p0, off_ao
      real(dp), intent(in) :: rad(:, :)  !! (chunk, 4) R, R', R'', R''' summed over primitives
      real(dp), intent(in) :: dd(:, :)   !! (chunk, 3) displacement from the shell's centre
      real(dp), intent(inout) :: ao(:, :)
      real(dp), intent(inout), optional :: grad(:, :, :)
      real(dp), intent(inout), optional :: hess(:, :, :)
      real(dp), intent(inout), optional :: deriv3(:, :, :)

      integer, parameter :: HESS_J(AO_HESS_COMP) = [1, 1, 1, 2, 2, 3]
      integer, parameter :: HESS_K(AO_HESS_COMP) = [1, 2, 3, 2, 3, 3]
      integer, parameter :: D3_I(AO_DERIV3_COMP) = [1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
      integer, parameter :: D3_J(AO_DERIV3_COMP) = [1, 1, 1, 2, 2, 3, 2, 2, 3, 3]
      integer, parameter :: D3_K(AO_DERIV3_COMP) = [1, 2, 3, 2, 3, 3, 2, 3, 3, 3]
      integer, parameter :: PAIR(3, 3) = reshape([1, 2, 3, 2, 4, 5, 3, 5, 6], [3, 3])

      integer :: s, col, p, ih, i, j, k
      real(dp) :: r0, r1, r2, r3, di, dj, dk, dij, dik, djk
      real(dp) :: ri, rj, rk, rij, rik, rjk, rijk

      do s = 1, n_fun
         col = off_ao + s
         do p = 1, np
            ao(p0 + p, col) = rad(p, 1)*a(p, s, 1)
         end do
         if (present(grad)) then
            do i = 1, 3
               do p = 1, np
                  grad(p0 + p, col, i) = a(p, s, 1 + i)*rad(p, 1) + a(p, s, 1)*rad(p, 2)*dd(p, i)
               end do
            end do
         end if
         if (present(hess)) then
            do ih = 1, AO_HESS_COMP
               i = HESS_J(ih)
               j = HESS_K(ih)
               dij = 0.0_dp
               if (i == j) dij = 1.0_dp
               do p = 1, np
                  hess(p0 + p, col, ih) = a(p, s, 4 + ih)*rad(p, 1) &
                                          + a(p, s, 1 + i)*rad(p, 2)*dd(p, j) &
                                          + a(p, s, 1 + j)*rad(p, 2)*dd(p, i) &
                                          + a(p, s, 1)*rad(p, 3)*dd(p, i)*dd(p, j) &
                                          + dij*a(p, s, 1)*rad(p, 2)
               end do
            end do
         end if
         if (present(deriv3)) then
            do ih = 1, AO_DERIV3_COMP
               i = D3_I(ih)
               j = D3_J(ih)
               k = D3_K(ih)
               dij = 0.0_dp
               dik = 0.0_dp
               djk = 0.0_dp
               if (i == j) dij = 1.0_dp
               if (i == k) dik = 1.0_dp
               if (j == k) djk = 1.0_dp
               do p = 1, np
                  r0 = rad(p, 1)
                  r1 = rad(p, 2)
                  r2 = rad(p, 3)
                  r3 = rad(p, 4)
                  di = dd(p, i)
                  dj = dd(p, j)
                  dk = dd(p, k)
                  ri = r1*di
                  rj = r1*dj
                  rk = r1*dk
                  rij = r2*di*dj + dij*r1
                  rik = r2*di*dk + dik*r1
                  rjk = r2*dj*dk + djk*r1
                  rijk = r3*di*dj*dk + r2*(dij*dk + dik*dj + djk*di)
                  deriv3(p0 + p, col, ih) = a(p, s, 10 + ih)*r0 &
                                            + a(p, s, 4 + PAIR(i, j))*rk &
                                            + a(p, s, 4 + PAIR(i, k))*rj &
                                            + a(p, s, 4 + PAIR(j, k))*ri &
                                            + a(p, s, 1 + i)*rjk &
                                            + a(p, s, 1 + j)*rik &
                                            + a(p, s, 1 + k)*rij &
                                            + a(p, s, 1)*rijk
               end do
            end do
         end if
      end do
   end subroutine store_shell

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
      !!
      !! `density` is whatever convention the caller's SCF uses -- the
      !! restricted path's D already carries two electrons per orbital, so this
      !! returns the total density for it and a spin density for an
      !! unrestricted D. **No factor is applied here.**
      real(dp), intent(in) :: ao(:, :)        !! (n_points, n_ao)
      real(dp), intent(in) :: density(:, :)   !! (n_ao, n_ao)
      real(dp), allocatable, intent(out) :: rho(:)
      real(dp), intent(in), optional :: ao_grad(:, :, :)
      real(dp), allocatable, intent(out), optional :: rho_grad(:, :)
         !! d rho / d r, as (n_points, 3). Both optionals or neither.
         !!
         !!     grad rho = 2 sum_uv D_uv chi_v grad chi_u
         !!
         !! The factor of two is the symmetry of D, not a spin factor.

      real(dp), allocatable, intent(out), optional :: tau(:)
         !! The kinetic energy density, for meta-GGA. Needs `ao_grad`.
         !!
         !!     tau(r) = 1/2 sum_d sum_uv D_uv d_d chi_u d_d chi_v
         !!
         !! The half is libxc's convention and PySCF's; the same quantity is
         !! defined without it elsewhere in the literature.

      real(dp), allocatable :: work(:, :), gwork(:, :)
      integer :: n_points, n_ao, ig, mu, id

      n_points = size(ao, 1)
      n_ao = size(ao, 2)

      allocate (work(n_points, n_ao), rho(n_points))
      call pic_gemm(ao, density, work)

      ! The row-wise dots walk the function index outermost, so every inner
      ! loop is one contiguous column: the other order strides by `n_points`
      ! on each step and was a quarter of a Kohn-Sham Hessian's samples. No
      ! threading here either -- every caller is already inside a parallel
      ! region over blocks of points, and a region nested there is a barrier
      ! that buys nothing.
      rho = 0.0_dp
      do mu = 1, n_ao
         do ig = 1, n_points
            rho(ig) = rho(ig) + work(ig, mu)*ao(ig, mu)
         end do
      end do

      if (present(tau)) then
         allocate (tau(n_points))
         tau = 0.0_dp
         if (.not. present(ao_grad)) then
            tau = huge(1.0_dp)
         else
            allocate (gwork(n_points, n_ao))
            do id = 1, 3
               call pic_gemm(ao_grad(:, :, id), density, gwork)
               do mu = 1, n_ao
                  do ig = 1, n_points
                     tau(ig) = tau(ig) + 0.5_dp*gwork(ig, mu)*ao_grad(ig, mu, id)
                  end do
               end do
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
            do id = 1, 3
               do mu = 1, n_ao
                  do ig = 1, n_points
                     rho_grad(ig, id) = rho_grad(ig, id) &
                                        + 2.0_dp*work(ig, mu)*ao_grad(ig, mu, id)
                  end do
               end do
            end do
         end if
      end if

      deallocate (work)
   end subroutine eval_rho

end module mqc_czt_ao
