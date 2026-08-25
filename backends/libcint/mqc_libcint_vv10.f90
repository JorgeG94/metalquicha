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
   implicit none
   private

   public :: vv10_nlc
   public :: VV10_RHO_THRESHOLD

   real(dp), parameter :: VV10_RHO_THRESHOLD = 1.0e-8_dp
      !! Points below this density are dropped from both grids.
      !!
      !! Not a tuning knob so much as the boundary of where the expression is
      !! defined: `omega_g` divides by rho squared and `kappa` takes rho to the
      !! sixth root, so vacuum points are simultaneously worthless and
      !! numerically hostile. PySCF uses the same 1e-8.

contains

   subroutine vv10_nlc(b, c, coords, rho, sigma, &
                       inner_coords, inner_rho, inner_sigma, inner_weights, &
                       exc, vrho, vsigma)
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

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: pi43, kvv, beta
      real(dp) :: w0, k_out, dw0_drho, dw0_dsigma, dk_drho, w0tmp
      real(dp) :: dx, dy, dz, r2, g, gp, gt, t, tt
      real(dp) :: f_sum, u_sum, w_sum
      real(dp) :: r_i, s_i
      real(dp), allocatable :: w0p(:), kp(:), rpw(:)
      integer :: n_out, n_in, i, j, n_kept
      integer, allocatable :: keep(:)

      n_out = size(rho)
      n_in = size(inner_rho)
      exc = 0.0_dp
      vrho = 0.0_dp
      vsigma = 0.0_dp

      pi43 = 4.0_dp*PI/3.0_dp
      kvv = b*1.5_dp*PI*(9.0_dp*PI)**(-1.0_dp/6.0_dp)
      beta = ((3.0_dp/(b*b))**0.75_dp)/32.0_dp

      !> The inner grid is precomputed once and then read n_out times, so the
      !> threshold pays for itself twice: it shrinks the array and it shortens
      !> every one of the outer loop's passes over it.
      allocate (keep(n_in))
      n_kept = 0
      do j = 1, n_in
         if (inner_rho(j) >= VV10_RHO_THRESHOLD) then
            n_kept = n_kept + 1
            keep(n_kept) = j
         end if
      end do
      if (n_kept == 0) return

      allocate (w0p(n_kept), kp(n_kept), rpw(n_kept))
      do j = 1, n_kept
         r_i = inner_rho(keep(j))
         s_i = inner_sigma(keep(j))
         w0tmp = c*(s_i/(r_i*r_i))**2
         w0p(j) = sqrt(w0tmp + pi43*r_i)
         kp(j) = kvv*r_i**(1.0_dp/6.0_dp)
         rpw(j) = r_i*inner_weights(keep(j))
      end do

      !> Every outer point is independent: it writes only its own three
      !> outputs and reads a shared inner grid that nothing mutates. The
      !> accumulators are per-point and so private. This loop is the whole cost
      !> of the term -- the two AO sweeps around it are linear in points and
      !> disappear against it -- so it is the only thing here worth threading.
      !$omp parallel do default(none) &
      !$omp    shared(n_out, n_kept, rho, sigma, coords, inner_coords, keep, &
      !$omp           w0p, kp, rpw, exc, vrho, vsigma, c, kvv, pi43, beta) &
      !$omp    private(i, j, r_i, s_i, w0, w0tmp, dw0_drho, dw0_dsigma, k_out, &
      !$omp            dk_drho, f_sum, u_sum, w_sum, dx, dy, dz, r2, g, gp, gt, t, tt) &
      !$omp    schedule(dynamic, 64)
      do i = 1, n_out
         if (rho(i) < VV10_RHO_THRESHOLD) cycle
         r_i = rho(i)
         s_i = sigma(i)

         w0tmp = c*(s_i/(r_i*r_i))**2
         w0 = sqrt(w0tmp + pi43*r_i)
         dw0_drho = (0.5_dp*pi43*r_i - 2.0_dp*w0tmp)/w0
         !> sigma appears in the denominator, and a point can have density with
         !> a vanishing gradient -- a nucleus, or any local extremum. w0tmp is
         !> proportional to sigma squared, so the limit is zero rather than
         !> singular; taking it explicitly avoids 0/0.
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
         do j = 1, n_kept
            dx = inner_coords(1, keep(j)) - coords(1, i)
            dy = inner_coords(2, keep(j)) - coords(2, i)
            dz = inner_coords(3, keep(j)) - coords(3, i)
            r2 = dx*dx + dy*dy + dz*dz
            gp = r2*w0p(j) + kp(j)
            g = r2*w0 + k_out
            gt = g + gp
            t = rpw(j)/(g*gp*gt)
            f_sum = f_sum + t
            tt = t*(1.0_dp/g + 1.0_dp/gt)
            u_sum = u_sum + tt
            w_sum = w_sum + tt*r2
         end do
         f_sum = -1.5_dp*f_sum

         exc(i) = beta + 0.5_dp*f_sum
         vrho(i) = beta + f_sum + 1.5_dp*(u_sum*dk_drho + w_sum*dw0_drho)
         vsigma(i) = 1.5_dp*w_sum*dw0_dsigma
      end do
      !$omp end parallel do
   end subroutine vv10_nlc

end module mqc_libcint_vv10
