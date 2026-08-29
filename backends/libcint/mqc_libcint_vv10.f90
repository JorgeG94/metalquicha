!! VV10 non-local correlation, as a double integral over the quadrature grid
module mqc_libcint_vv10
   !! Vydrov and Van Voorhis, J. Chem. Phys. 133, 244103 (2010).
   !!
   !!     E_c^nl = int rho(r) [ 1/2 int rho(r') Phi(r,r') dr' + beta ] dr
   !!
   !! **This is not a functional of the density at a point**, which is why it
   !! cannot ride along inside libxc the way everything else here does. libxc
   !! supplies the semilocal half of a `-V` functional and hands back the two
   !! parameters `b` and `c`; the double integral is the program's job. Skipping
   !! it does not approximate the functional -- on water/STO-3G it moves
   !! wB97X-V by 43 mHa, 27 kcal/mol -- so `mqc_libcint_xc` refuses a `-V`
   !! functional outright unless this is available.
   !!
   !! **Cost is the design constraint.** The pair sum is O(N_out * N_in) in
   !! grid points, per SCF iteration, and mqc's default grid puts 33,704 points
   !! on water and 84,712 on benzene. Evaluating this on the exchange grid is
   !! not viable, so the caller is expected to pass a coarser one for the inner
   !! sum, and the density threshold below removes most of what is left: a
   !! molecular grid is mostly vacuum, and a point with no density contributes
   !! nothing to either side of the integral.
   !!
   !! The arithmetic follows PySCF's `_vv10nlc` term for term, deliberately, so
   !! that a disagreement is a porting bug rather than a difference of
   !! convention. Its published pure-Python inner loop is the reference this was
   !! written from and the unit test checks against.
   use pic_types, only: dp
   use mqc_physical_constants, only: PI
   implicit none
   private

   public :: vv10_nlc
   public :: vv10_hessian_kernel
   public :: VV10_RHO_THRESHOLD

   real(dp), parameter :: VV10_RHO_THRESHOLD = 1.0e-8_dp

   !> Beyond this separation, in Bohr, an inner point is not summed. Negative
   !> disables the cutoff and restores the full double sum.
   !>
   !> The kernel falls off as r^-6: at large separation `g -> r^2 w0`,
   !> `gp -> r^2 w0p` and `gt -> r^2 (w0 + w0p)`, so `t = rpw/(g gp gt)` goes as
   !> `rpw / r^6`. The number of points at a separation grows as r^2, so the
   !> tail this drops integrates as r^-3 -- fast, but not so fast that the
   !> radius can be picked by eye. It is measured, not assumed.
   real(dp), parameter, public :: VV10_CUTOFF = -1.0_dp
      !! Points below this density are dropped from both grids.
      !!
      !! Not a tuning knob so much as the boundary of where the expression is
      !! defined: `omega_g` divides by rho squared and `kappa` takes rho to the
      !! sixth root, so vacuum points are simultaneously worthless and
      !! numerically hostile. PySCF uses the same 1e-8.

contains

   subroutine vv10_nlc(b, c, coords, rho, sigma, &
                       inner_coords, inner_rho, inner_sigma, inner_weights, &
                       exc, vrho, vsigma, dedw, fexp, &
                       hess_u, hess_w, hess_a, hess_b, hess_c, hess_e, &
                       domega_drho, domega_dgamma, d2omega_drho2, &
                       d2omega_dgamma2, d2omega_drho_dgamma, &
                       dkappa_drho, d2kappa_drho2)
      !! VV10's energy density and potential at each outer grid point
      !!
      !! `exc` is the energy *per electron*, so the caller forms the energy the
      !! same way it does for every other functional:
      !!
      !!     E_nl = sum_g w_g rho_g exc_g
      !!
      !! and `vrho`, `vsigma` add into the quantities `xc_grid_gga_quantities`
      !! already returns. That is the whole point of this signature: VV10
      !! becomes another contribution to numbers the Fock build already
      !! consumes, rather than a second path through the SCF.
      real(dp), intent(in) :: b, c            !! libxc's nlc_b and nlc_c
      real(dp), intent(in) :: coords(:, :)    !! (3, n_out), Bohr
      real(dp), intent(in) :: rho(:)          !! (n_out)
      real(dp), intent(in) :: sigma(:)        !! (n_out), |grad rho|^2
      real(dp), intent(in) :: inner_coords(:, :)   !! (3, n_in)
      real(dp), intent(in) :: inner_rho(:)         !! (n_in)
      real(dp), intent(in) :: inner_sigma(:)       !! (n_in)
      real(dp), intent(in) :: inner_weights(:)     !! (n_in)
      real(dp), intent(out) :: exc(:)         !! (n_out)
      real(dp), intent(out) :: vrho(:)        !! (n_out), dE/drho
      real(dp), intent(out) :: vsigma(:)      !! (n_out), dE/dsigma
      real(dp), intent(out), optional :: dedw(:)
         !! (n_out), dE/dw_i -- the derivative with respect to the *weight* of
         !! a point, which a nuclear gradient needs because the Becke partition
         !! moves with the nuclei.
         !!
         !! **It is not `rho*exc`.** For a semilocal functional it would be:
         !! the energy is a sum of independent per-point terms and the weight
         !! of one multiplies one of them. Here point `k` appears twice, once
         !! as the outer point and once inside every other point's inner sum,
         !! and the kernel's symmetry makes those two contributions equal. So
         !! this carries `beta + f` where `exc` carries `beta + f/2` -- the
         !! same doubling that distinguishes `vrho` from `exc` just above.
      real(dp), intent(out), optional :: fexp(:, :)
         !! (3, n_out), dE/dr_i at fixed density and weights, per unit weight
         !!
         !! The energy depends on *where the points are* and not only on what
         !! the density is there: `r^2` sits inside the kernel. Nothing in a
         !! semilocal functional does this, which is why a gradient assembled
         !! from the usual three terms is incomplete here by about 4e-4 on
         !! water and does not converge away with the grid.
         !!
         !! Per unit weight, so the caller multiplies by `w_i` exactly as it
         !! does for `exc`, and sums onto the atom that owns the point.
      real(dp), intent(out), optional :: hess_u(:), hess_w(:)
         !! (n_out) each. The pair sums a second derivative needs, in PySCF's
         !! `VXC_vv10nlc_hessian_eval_UWABCE` normalisation, so that every
         !! formula in the Hessian papers can be transcribed rather than
         !! re-derived. With Phi the pair kernel and g, g', g_s = g + g' its
         !! denominators evaluated at the outer point's omega and kappa:
         !!
         !!     U = -sum_j w_j rho_j Phi (1/g + 1/g_s)
         !!     W = -sum_j w_j rho_j Phi (1/g + 1/g_s) r^2
         !!
         !! These are the same sums `vrho` is built from, scaled: PySCF's
         !! U is `1.5*u_sum` here, W is `1.5*w_sum`.
      real(dp), intent(out), optional :: hess_a(:), hess_b(:), hess_c(:)
         !! (n_out) each. The genuinely new pair sums, quadratic in the
         !! denominators -- what differentiating U and W a second time
         !! produces:
         !!
         !!     A = 2 sum_j w_j rho_j Phi (1/g^2 + 1/(g g_s) + 1/g_s^2)
         !!     B = A's sum with one factor of r^2
         !!     C = A's sum with two
      real(dp), intent(out), optional :: hess_e(:)
         !! (n_out). E = sum_j w_j rho_j Phi -- the inner sum itself, which is
         !! `2*(exc - beta)`. Returned anyway because PySCF's Hessian formulas
         !! consume it under this name and normalisation.
      real(dp), intent(out), optional :: domega_drho(:), domega_dgamma(:)
      real(dp), intent(out), optional :: d2omega_drho2(:), d2omega_dgamma2(:)
      real(dp), intent(out), optional :: d2omega_drho_dgamma(:)
      real(dp), intent(out), optional :: dkappa_drho(:), d2kappa_drho2(:)
         !! (n_out) each: analytic derivatives of omega_0 and kappa at the
         !! outer points, with respect to rho and gamma = sigma themselves --
         !! **not** premultiplied by rho the way the inline `dw0_drho` and
         !! `dw0_dsigma` below are. PySCF's convention again, and the test
         !! pins the relation: `rho*domega_drho` must rebuild `dw0_drho`.

      real(dp) :: pi43, kvv, beta
      real(dp) :: w0, k_out, dw0_drho, dw0_dsigma, dk_drho, w0tmp
      real(dp) :: dx, dy, dz, r2, g, gp, gt, t, tt, q
      real(dp) :: fx, fy, fz
      logical :: want_grad, want_hess, want_curv
      real(dp) :: f_sum, u_sum, w_sum
      real(dp) :: ta, a_sum, b_sum, c_sum
      real(dp) :: rho_1, rho_3, rho_4, rho_5, gamma2, omega_1, omega_3
      real(dp) :: r_i, s_i
      real(dp), allocatable :: w0p(:), kp(:), rpw(:)
      real(dp), allocatable :: cx(:), cy(:), cz(:)
      real(dp) :: r2_cut
      integer :: n_out, n_in, i, j, n_kept
      integer, allocatable :: keep(:)

      n_out = size(rho)
      n_in = size(inner_rho)
      exc = 0.0_dp
      vrho = 0.0_dp
      vsigma = 0.0_dp
      want_grad = present(dedw) .or. present(fexp)
      if (present(dedw)) dedw = 0.0_dp
      if (present(fexp)) fexp = 0.0_dp
      ! Guarded separately from `want_grad` because these are Hessian-only and
      ! the potential path runs every SCF iteration: A, B and C add work to
      ! the inner loop, which is the whole cost of this routine.
      want_hess = present(hess_a) .or. present(hess_b) .or. present(hess_c)
      want_curv = present(domega_drho) .or. present(domega_dgamma) .or. &
                  present(d2omega_drho2) .or. present(d2omega_dgamma2) .or. &
                  present(d2omega_drho_dgamma) .or. present(dkappa_drho) .or. &
                  present(d2kappa_drho2)
      if (present(hess_u)) hess_u = 0.0_dp
      if (present(hess_w)) hess_w = 0.0_dp
      if (present(hess_a)) hess_a = 0.0_dp
      if (present(hess_b)) hess_b = 0.0_dp
      if (present(hess_c)) hess_c = 0.0_dp
      if (present(hess_e)) hess_e = 0.0_dp
      if (present(domega_drho)) domega_drho = 0.0_dp
      if (present(domega_dgamma)) domega_dgamma = 0.0_dp
      if (present(d2omega_drho2)) d2omega_drho2 = 0.0_dp
      if (present(d2omega_dgamma2)) d2omega_dgamma2 = 0.0_dp
      if (present(d2omega_drho_dgamma)) d2omega_drho_dgamma = 0.0_dp
      if (present(dkappa_drho)) dkappa_drho = 0.0_dp
      if (present(d2kappa_drho2)) d2kappa_drho2 = 0.0_dp

      pi43 = 4.0_dp*PI/3.0_dp
      kvv = b*1.5_dp*PI*(9.0_dp*PI)**(-1.0_dp/6.0_dp)
      beta = ((3.0_dp/(b*b))**0.75_dp)/32.0_dp

      ! The inner grid is precomputed once and then read n_out times, so the
      ! threshold pays for itself twice: it shrinks the array and it shortens
      ! every one of the outer loop's passes over it.
      allocate (keep(n_in))
      n_kept = 0
      do j = 1, n_in
         if (inner_rho(j) >= VV10_RHO_THRESHOLD) then
            n_kept = n_kept + 1
            keep(n_kept) = j
         end if
      end do
      if (n_kept == 0) return

      ! The coordinates are gathered alongside, rather than reached through
      ! `keep` in the inner loop. That loop is the whole cost of VV10, and an
      ! indirection inside it is an extra load per pair and a barrier to
      ! vectorising the three differences.
      allocate (w0p(n_kept), kp(n_kept), rpw(n_kept))
      allocate (cx(n_kept), cy(n_kept), cz(n_kept))
      do j = 1, n_kept
         r_i = inner_rho(keep(j))
         s_i = inner_sigma(keep(j))
         w0tmp = c*(s_i/(r_i*r_i))**2
         w0p(j) = sqrt(w0tmp + pi43*r_i)
         kp(j) = kvv*r_i**(1.0_dp/6.0_dp)
         rpw(j) = r_i*inner_weights(keep(j))
         cx(j) = inner_coords(1, keep(j))
         cy(j) = inner_coords(2, keep(j))
         cz(j) = inner_coords(3, keep(j))
      end do

      r2_cut = huge(1.0_dp)
      if (VV10_CUTOFF > 0.0_dp) r2_cut = VV10_CUTOFF*VV10_CUTOFF

      ! Every outer point is independent: it writes only its own three
      ! outputs and reads a shared inner grid that nothing mutates. The
      ! accumulators are per-point and so private. This loop is the whole cost
      ! of the term -- the two AO sweeps around it are linear in points and
      ! disappear against it -- so it is the only thing here worth threading.
      !$omp parallel do default(none) &
      !$omp    shared(n_out, n_kept, rho, sigma, coords, inner_coords, keep, &
      !$omp           w0p, kp, rpw, cx, cy, cz, r2_cut, exc, vrho, vsigma, c, kvv, pi43, beta, &
      !$omp           want_grad, dedw, fexp, want_hess, want_curv, &
      !$omp           hess_u, hess_w, hess_a, hess_b, hess_c, hess_e, &
      !$omp           domega_drho, domega_dgamma, d2omega_drho2, &
      !$omp           d2omega_dgamma2, d2omega_drho_dgamma, &
      !$omp           dkappa_drho, d2kappa_drho2) &
      !$omp    private(i, j, r_i, s_i, w0, w0tmp, dw0_drho, dw0_dsigma, k_out, &
      !$omp            dk_drho, f_sum, u_sum, w_sum, dx, dy, dz, r2, g, gp, gt, t, tt, &
      !$omp            q, fx, fy, fz, ta, a_sum, b_sum, c_sum, &
      !$omp            rho_1, rho_3, rho_4, rho_5, gamma2, omega_1, omega_3) &
      !$omp    schedule(dynamic, 64)
      do i = 1, n_out
         if (rho(i) < VV10_RHO_THRESHOLD) cycle
         r_i = rho(i)
         s_i = sigma(i)

         w0tmp = c*(s_i/(r_i*r_i))**2
         w0 = sqrt(w0tmp + pi43*r_i)
         dw0_drho = (0.5_dp*pi43*r_i - 2.0_dp*w0tmp)/w0
         ! sigma appears in the denominator, and a point can have density with
         ! a vanishing gradient -- a nucleus, or any local extremum. w0tmp is
         ! proportional to sigma squared, so the limit is zero rather than
         ! singular; taking it explicitly avoids 0/0.
         if (s_i > 0.0_dp) then
            dw0_dsigma = w0tmp*r_i/(s_i*w0)
         else
            dw0_dsigma = 0.0_dp
         end if
         k_out = kvv*r_i**(1.0_dp/6.0_dp)
         dk_drho = k_out/6.0_dp

         f_sum = 0.0_dp
         u_sum = 0.0_dp
         w_sum = 0.0_dp
         a_sum = 0.0_dp
         b_sum = 0.0_dp
         c_sum = 0.0_dp
         fx = 0.0_dp
         fy = 0.0_dp
         fz = 0.0_dp
         do j = 1, n_kept
            dx = cx(j) - coords(1, i)
            dy = cy(j) - coords(2, i)
            dz = cz(j) - coords(3, i)
            r2 = dx*dx + dy*dy + dz*dz
            ! Three divisions and a dozen multiplies are skipped for the price
            ! of a compare. Visiting the pair at all is still O(n^2); making it
            ! sub-quadratic needs the inner list binned in space, which this is
            ! the prerequisite for rather than a substitute for.
            if (r2 > r2_cut) cycle
            gp = r2*w0p(j) + kp(j)
            g = r2*w0 + k_out
            gt = g + gp
            t = rpw(j)/(g*gp*gt)
            f_sum = f_sum + t
            tt = t*(1.0_dp/g + 1.0_dp/gt)
            u_sum = u_sum + tt
            w_sum = w_sum + tt*r2
            ! Differentiating 1/g and 1/(g+g') once more, hence the squares
            ! and the cross term. Only the second derivative reads these.
            if (want_hess) then
               ta = t*(1.0_dp/(g*g) + 1.0_dp/(g*gt) + 1.0_dp/(gt*gt))
               a_sum = a_sum + ta
               b_sum = b_sum + ta*r2
               c_sum = c_sum + ta*r2*r2
            end if
            ! The kernel's own dependence on the separation. `g` and `gp` each
            ! carry `r^2`, weighted by that end's `w0`, and `gt` carries both.
            if (want_grad) then
               q = t*(w0/g + w0p(j)/gp + (w0 + w0p(j))/gt)
               fx = fx + q*dx
               fy = fy + q*dy
               fz = fz + q*dz
            end if
         end do
         f_sum = -1.5_dp*f_sum

         exc(i) = beta + 0.5_dp*f_sum
         vrho(i) = beta + f_sum + 1.5_dp*(u_sum*dk_drho + w_sum*dw0_drho)
         vsigma(i) = 1.5_dp*w_sum*dw0_dsigma
         if (present(dedw)) dedw(i) = r_i*(beta + f_sum)
         if (present(fexp)) then
            ! d/dr_i of the pair sum: two from differentiating a squared
            ! separation, and another from the outer point appearing in both
            ! roles, against `dx` which is measured inner minus outer.
            fexp(1, i) = -3.0_dp*r_i*fx
            fexp(2, i) = -3.0_dp*r_i*fy
            fexp(3, i) = -3.0_dp*r_i*fz
         end if
         ! PySCF's signs and scales, term for term: its per-pair kernel is
         ! `-1.5*t`, its U and W carry an overall minus, its A, B, C a factor
         ! of two. Written out so that a disagreement with its Hessian is a
         ! porting bug and not a convention.
         if (present(hess_u)) hess_u(i) = 1.5_dp*u_sum
         if (present(hess_w)) hess_w(i) = 1.5_dp*w_sum
         if (present(hess_a)) hess_a(i) = -3.0_dp*a_sum
         if (present(hess_b)) hess_b(i) = -3.0_dp*b_sum
         if (present(hess_c)) hess_c(i) = -3.0_dp*c_sum
         if (present(hess_e)) hess_e(i) = f_sum
         if (want_curv) then
            ! Derivatives of omega_0 and kappa with respect to rho and
            ! gamma = sigma, PySCF's `_omega_derivative` verbatim. Every
            ! expression is polynomial in gamma -- gamma multiplies, it never
            ! divides -- so the sigma -> 0 limit that `dw0_dsigma` above has
            ! to take explicitly is automatic here. Writing these via
            ! `w0tmp/sigma` instead would reintroduce that 0/0 at every
            ! nucleus and density extremum.
            rho_1 = 1.0_dp/r_i
            rho_3 = rho_1*rho_1*rho_1
            rho_4 = rho_3*rho_1
            rho_5 = rho_4*rho_1
            gamma2 = s_i*s_i
            omega_1 = 1.0_dp/w0
            omega_3 = omega_1/(w0*w0)
            if (present(domega_drho)) domega_drho(i) = &
               0.5_dp*(pi43 - 4.0_dp*c*gamma2*rho_5)*omega_1
            if (present(domega_dgamma)) domega_dgamma(i) = &
               c*s_i*rho_4*omega_1
            if (present(d2omega_drho2)) d2omega_drho2(i) = &
               (-0.25_dp*pi43*pi43 + 12.0_dp*pi43*c*gamma2*rho_5 &
                + 6.0_dp*c*c*gamma2*gamma2*rho_5*rho_5)*omega_3
            if (present(d2omega_dgamma2)) d2omega_dgamma2(i) = &
               pi43*c*rho_3*omega_3
            if (present(d2omega_drho_dgamma)) d2omega_drho_dgamma(i) = &
               -c*s_i*(4.5_dp*pi43*rho_4 + 2.0_dp*c*gamma2*rho_4*rho_5)*omega_3
            if (present(dkappa_drho)) dkappa_drho(i) = k_out/(6.0_dp*r_i)
            if (present(d2kappa_drho2)) d2kappa_drho2(i) = &
               -5.0_dp*k_out/(36.0_dp*r_i*r_i)
         end if
      end do
      !$omp end parallel do
   end subroutine vv10_nlc

   subroutine vv10_hessian_kernel(b, c, coords, rho, sigma, weights, &
                                  hess_u, hess_w, hess_a, hess_b, hess_c, &
                                  domega_drho, domega_dgamma, dkappa_drho, &
                                  d2omega_drho2, d2omega_dgamma2, &
                                  d2omega_drho_dgamma, d2kappa_drho2, &
                                  rho_t, gamma_t, f_rho_t, f_gamma_t)
      !! The VV10 kernel applied to a batch of density perturbations
      !!
      !! PySCF's `VXC_vv10nlc_hessian_eval_f_t`, transcribed term for term.
      !! For each trial pair `(rho_t, gamma_t)` -- in a Hessian these are
      !! `drho/dA` and `dgamma/dA` for one nuclear perturbation, but nothing
      !! here knows that -- it returns the second functional derivative of
      !! `E_nl` applied to it:
      !!
      !!     f_rho_t(i)   = sum_j w_j [ f_rr(i,j) rho_t(j) + f_rg(i,j) gamma_t(j) ]
      !!                  + f_rr_ii rho_t(i) + f_rg_ii gamma_t(i)
      !!
      !! and its gamma-channel twin. The double sum is the kernel's genuinely
      !! non-local part -- the potential at `i` feels the density at `j`
      !! through the pair kernel Phi -- and the unpaired diagonal is the
      !! potential's dependence on its *own* point's rho and gamma through
      !! omega and kappa, which is where the pair sums U..C and the second
      !! derivatives of omega and kappa enter. The inner quadrature weight
      !! `w_j` is applied here; the outer `w_i` is the caller's, exactly as
      !! PySCF weights `f_rho_t` only when contracting it.
      !!
      !! One grid plays both roles, unlike `vv10_nlc`'s two: the Hessian's
      !! trial densities live on the grid the intermediates were built on,
      !! and passing them separately would only invite a mismatch.
      real(dp), intent(in) :: b, c              !! libxc's nlc_b and nlc_c
      real(dp), intent(in) :: coords(:, :)      !! (3, n), Bohr
      real(dp), intent(in) :: rho(:), sigma(:), weights(:)
      real(dp), intent(in) :: hess_u(:), hess_w(:)
      real(dp), intent(in) :: hess_a(:), hess_b(:), hess_c(:)
         !! `vv10_nlc`'s pair sums, in the same normalisation it returns them
      real(dp), intent(in) :: domega_drho(:), domega_dgamma(:), dkappa_drho(:)
      real(dp), intent(in) :: d2omega_drho2(:), d2omega_dgamma2(:)
      real(dp), intent(in) :: d2omega_drho_dgamma(:), d2kappa_drho2(:)
      real(dp), intent(in) :: rho_t(:, :), gamma_t(:, :)      !! (n_trial, n)
      real(dp), intent(out) :: f_rho_t(:, :), f_gamma_t(:, :)  !! (n_trial, n)

      real(dp) :: pi43, kvv
      real(dp) :: r_i, s_i, w0tmp, w0_i, k_i, dor_i, dog_i, dkr_i
      real(dp) :: dx, dy, dz, r2, g, gp, gt, g1, gp1, gt1, phi
      real(dp) :: rdg_i, rdg_j, d2phi, f_rr, f_rg, f_gr, f_gg
      integer :: n, n_trial, n_kept, i, ii, jj, it
      integer, allocatable :: keep(:)
      real(dp), allocatable :: w0p(:), kp(:)

      n = size(rho)
      n_trial = size(rho_t, 1)
      f_rho_t = 0.0_dp
      f_gamma_t = 0.0_dp

      pi43 = 4.0_dp*PI/3.0_dp
      kvv = b*1.5_dp*PI*(9.0_dp*PI)**(-1.0_dp/6.0_dp)

      ! The same threshold, on both roles at once: PySCF removes these points
      ! from the quadrature before its kernel ever runs, and omega and kappa
      ! are not defined below it anyway.
      allocate (keep(n))
      n_kept = 0
      do i = 1, n
         if (rho(i) >= VV10_RHO_THRESHOLD) then
            n_kept = n_kept + 1
            keep(n_kept) = i
         end if
      end do
      if (n_kept == 0) return

      allocate (w0p(n_kept), kp(n_kept))
      do jj = 1, n_kept
         r_i = rho(keep(jj))
         s_i = sigma(keep(jj))
         w0tmp = c*(s_i/(r_i*r_i))**2
         w0p(jj) = sqrt(w0tmp + pi43*r_i)
         kp(jj) = kvv*r_i**(1.0_dp/6.0_dp)
      end do

      ! Outer points are independent, as in `vv10_nlc`; each iteration writes
      ! only its own two output columns, so they can be accumulated in place.
      !$omp parallel do default(none) &
      !$omp    shared(n_kept, n_trial, keep, coords, rho, weights, w0p, kp, &
      !$omp           hess_u, hess_w, hess_a, hess_b, hess_c, &
      !$omp           domega_drho, domega_dgamma, dkappa_drho, &
      !$omp           d2omega_drho2, d2omega_dgamma2, d2omega_drho_dgamma, &
      !$omp           d2kappa_drho2, rho_t, gamma_t, f_rho_t, f_gamma_t) &
      !$omp    private(i, ii, jj, it, r_i, w0_i, k_i, dor_i, dog_i, dkr_i, &
      !$omp            dx, dy, dz, r2, g, gp, gt, g1, gp1, gt1, phi, &
      !$omp            rdg_i, rdg_j, d2phi, f_rr, f_rg, f_gr, f_gg) &
      !$omp    schedule(dynamic, 16)
      do ii = 1, n_kept
         i = keep(ii)
         r_i = rho(i)
         w0_i = w0p(ii)
         k_i = kp(ii)
         dor_i = domega_drho(i)
         dog_i = domega_dgamma(i)
         dkr_i = dkappa_drho(i)

         ! The pair part, `j == i` included: that term is the two points
         ! coinciding inside the double integral, at `r2 = 0`, and is distinct
         ! from the diagonal added after the loop.
         do jj = 1, n_kept
            dx = coords(1, keep(jj)) - coords(1, i)
            dy = coords(2, keep(jj)) - coords(2, i)
            dz = coords(3, keep(jj)) - coords(3, i)
            r2 = dx*dx + dy*dy + dz*dz
            g = r2*w0_i + k_i
            gp = r2*w0p(jj) + kp(jj)
            gt = g + gp
            g1 = 1.0_dp/g
            gp1 = 1.0_dp/gp
            gt1 = 1.0_dp/gt
            phi = -1.5_dp*g1*gp1*gt1

            ! rho d(g)/drho at each end, and the second derivative of Phi with
            ! respect to its two denominators, over Phi.
            rdg_i = r_i*(r2*dor_i + dkr_i)
            rdg_j = rho(keep(jj))*(r2*domega_drho(keep(jj)) + dkappa_drho(keep(jj)))
            d2phi = 2.0_dp*(gt1*gt1 + g1*gp1)

            f_rr = phi*(rdg_i*rdg_j*d2phi - rdg_i*(gt1 + g1) &
                        - rdg_j*(gt1 + gp1) + 1.0_dp)
            f_gr = r_i*dog_i*r2*phi*(rdg_j*d2phi - (gt1 + g1))
            f_rg = rho(keep(jj))*domega_dgamma(keep(jj))*r2*phi &
                   *(rdg_i*d2phi - (gt1 + gp1))
            f_gg = r_i*rho(keep(jj))*dog_i*domega_dgamma(keep(jj)) &
                   *r2*r2*phi*d2phi

            do it = 1, n_trial
               f_rho_t(it, i) = f_rho_t(it, i) + weights(keep(jj)) &
                                *(f_rr*rho_t(it, keep(jj)) + f_rg*gamma_t(it, keep(jj)))
               f_gamma_t(it, i) = f_gamma_t(it, i) + weights(keep(jj)) &
                                  *(f_gr*rho_t(it, keep(jj)) + f_gg*gamma_t(it, keep(jj)))
            end do
         end do

         ! The diagonal: omega and kappa at `i` moving with rho and gamma
         ! there, felt through every pair `i` participates in -- which is
         ! exactly what the pair sums U..C already integrated, `w_j` and all,
         ! so no quadrature weight is applied here.
         f_rr = 2.0_dp*dor_i*hess_w(i) + 2.0_dp*dkr_i*hess_u(i) &
                + r_i*(d2omega_drho2(i)*hess_w(i) + d2kappa_drho2(i)*hess_u(i) &
                       + dkr_i*dkr_i*hess_a(i) + dor_i*dor_i*hess_c(i) &
                       + 2.0_dp*dor_i*dkr_i*hess_b(i))
         f_gr = dog_i*hess_w(i) + r_i*(d2omega_drho_dgamma(i)*hess_w(i) &
                                       + dog_i*(dkr_i*hess_b(i) + dor_i*hess_c(i)))
         f_rg = f_gr
         f_gg = r_i*(d2omega_dgamma2(i)*hess_w(i) + dog_i*dog_i*hess_c(i))

         do it = 1, n_trial
            f_rho_t(it, i) = f_rho_t(it, i) &
                             + f_rr*rho_t(it, i) + f_rg*gamma_t(it, i)
            f_gamma_t(it, i) = f_gamma_t(it, i) &
                               + f_gr*rho_t(it, i) + f_gg*gamma_t(it, i)
         end do
      end do
      !$omp end parallel do
   end subroutine vv10_hessian_kernel

end module mqc_libcint_vv10
