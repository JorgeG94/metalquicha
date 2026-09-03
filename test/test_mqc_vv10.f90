!! The VV10 non-local correlation kernel, against its reference implementation
!!
!! The numbers below were produced by PySCF's own `_vv10nlc` arithmetic, run on
!! the five points declared here. That is a deliberate choice of reference: the
!! risk in this kernel is a porting mistake -- a sign, a power, a term dropped
!! from the chain rule -- and comparing against the source it was ported from
!! isolates exactly that. Comparing against a physical energy instead would mix
!! a porting bug with a grid difference and diagnose neither.
!!
!! The five points are deliberately asymmetric, with no repeated densities,
!! gradients or spacings. VV10's kernel is symmetric in its two arguments, so a
!! symmetric input can hide an index swapped between the inner and outer grids.
!!
!! `sigma` at point 3 is small but non-zero; the zero-gradient limit in
!! `dw0_dsigma` is exercised separately, because it is a 0/0 that a physical
!! grid reaches at every nucleus.
module test_mqc_vv10
   use pic_types, only: dp
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mqc_error, only: error_t
   use mqc_czt_vv10, only: vv10_nlc, VV10_RHO_THRESHOLD
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available, &
                         ensure_nlc_grid
   use mqc_czt_ao, only: eval_ao_block, eval_rho
   implicit none
   private

   public :: collect_vv10

   integer, parameter :: NPT = 5

   !! wB97X-V's parameters, as libxc reports them through xc_f03_nlc_coef
   real(dp), parameter :: B_VV = 5.9_dp
   real(dp), parameter :: C_VV = 0.0093_dp

   real(dp), parameter :: COORDS(3, NPT) = reshape([ &
                                                   0.0_dp, 0.0_dp, 0.0_dp, &
                                                   1.2_dp, -0.4_dp, 0.7_dp, &
                                                   -0.9_dp, 1.1_dp, -1.3_dp, &
                                                   0.3_dp, 0.3_dp, 2.0_dp, &
                                                   -2.0_dp, -1.5_dp, 0.2_dp], [3, NPT])
   real(dp), parameter :: RHO(NPT) = [0.30_dp, 0.12_dp, 0.55_dp, 0.04_dp, 0.21_dp]
   real(dp), parameter :: SIGMA(NPT) = [0.05_dp, 0.31_dp, 0.02_dp, 0.44_dp, 0.10_dp]
   real(dp), parameter :: WEIGHTS(NPT) = [0.70_dp, 1.30_dp, 0.90_dp, 1.10_dp, 0.60_dp]

   real(dp), parameter :: EXC_REF(NPT) = [ &
                          4.88389639834253e-03_dp, 4.90063660795110e-03_dp, &
                          4.89201157661231e-03_dp, 4.94888481470596e-03_dp, &
                          4.92099962904201e-03_dp]
   real(dp), parameter :: VRHO_REF(NPT) = [ &
                          4.85153333869646e-03_dp, 4.80171426879059e-03_dp, &
                          4.85880635020367e-03_dp, 4.92824897791482e-03_dp, &
                          4.90135859343933e-03_dp]
   real(dp), parameter :: VSIGMA_REF(NPT) = [ &
                          4.74823422798564e-07_dp, 1.16540461258463e-05_dp, &
                          8.59841311166584e-09_dp, 4.16019979984147e-07_dp, &
                          2.97672128300922e-06_dp]
   real(dp), parameter :: E_NL_REF = 5.04946017002175e-03_dp

   !! The reference carries fourteen digits, so this is close to what double
   !! precision can express over a sum of this length. A looser threshold would
   !! accept a dropped term small enough to matter later.
   real(dp), parameter :: THR = 1.0e-12_dp

contains

   subroutine collect_vv10(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("vv10 energy density matches the reference", test_exc), &
                  new_unittest("vv10 potential matches the reference", test_potential), &
                  new_unittest("vv10 integrates to the reference energy", test_energy), &
                  new_unittest("vv10 survives a vanishing gradient", test_zero_gradient), &
                  new_unittest("vv10 ignores points below threshold", test_threshold), &
                  new_unittest("vv10 hessian intermediates rebuild the potential", &
                               test_rebuild_synthetic), &
                  new_unittest("vv10 intermediates rebuild the potential on a real grid", &
                               test_rebuild_real_grid) &
                  ]
   end subroutine collect_vv10

   subroutine rebuild_check(b, c, coords, rho, sigma, weights, thr, thr_rel, what, error)
      !! The identity that pins the intermediates' conventions
      !!
      !! PySCF's Hessian never touches `vrho` and `vsigma`; it rebuilds the
      !! potential from the intermediates as
      !!
      !!     f_rho   = beta + E + rho*(dkappa_drho*U + domega_drho*W)
      !!     f_gamma = rho*domega_dgamma*W
      !!
      !! and every later formula assumes those are the same numbers. So the
      !! rebuild must reproduce the validated `vrho` and `vsigma` to rounding:
      !! any normalisation or sign slip in U, W, E or the omega and kappa
      !! derivatives lands here as a finite disagreement, before anything
      !! downstream can inherit it.
      real(dp), intent(in) :: b, c
      real(dp), intent(in) :: coords(:, :), rho(:), sigma(:), weights(:)
      real(dp), intent(in) :: thr, thr_rel
      character(len=*), intent(in) :: what
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: beta, worst_rho, worst_gamma
      real(dp), allocatable :: exc(:), vrho(:), vsigma(:)
      real(dp), allocatable :: u(:), w(:), a(:), bb(:), cc(:), e(:)
      real(dp), allocatable :: dodr(:), dodg(:), d2odr2(:), d2odg2(:), d2odrdg(:)
      real(dp), allocatable :: dkdr(:), d2kdr2(:)
      real(dp) :: f_rho, f_gamma
      integer :: n, i, n_live

      n = size(rho)
      allocate (exc(n), vrho(n), vsigma(n))
      allocate (u(n), w(n), a(n), bb(n), cc(n), e(n))
      allocate (dodr(n), dodg(n), d2odr2(n), d2odg2(n), d2odrdg(n))
      allocate (dkdr(n), d2kdr2(n))

      call vv10_nlc(b, c, coords, rho, sigma, &
                    coords, rho, sigma, weights, exc, vrho, vsigma, &
                    hess_u=u, hess_w=w, hess_a=a, hess_b=bb, hess_c=cc, &
                    hess_e=e, domega_drho=dodr, domega_dgamma=dodg, &
                    d2omega_drho2=d2odr2, d2omega_dgamma2=d2odg2, &
                    d2omega_drho_dgamma=d2odrdg, dkappa_drho=dkdr, &
                    d2kappa_drho2=d2kdr2)

      beta = ((3.0_dp/(b*b))**0.75_dp)/32.0_dp
      worst_rho = 0.0_dp
      worst_gamma = 0.0_dp
      n_live = 0
      do i = 1, n
         ! Below the density threshold everything is zeroed, including vrho,
         ! while the rebuild would still carry beta; the identity only holds
         ! where the kernel is evaluated.
         if (rho(i) < VV10_RHO_THRESHOLD) cycle
         n_live = n_live + 1
         f_rho = beta + e(i) + rho(i)*(dkdr(i)*u(i) + dodr(i)*w(i))
         f_gamma = rho(i)*dodg(i)*w(i)
         worst_rho = max(worst_rho, abs(f_rho - vrho(i)))
         ! Relatively, as `test_potential` compares vsigma: near the density
         ! threshold vsigma itself grows like rho**(-3), so a fixed absolute
         ! band would be tight at a bond midpoint and vacuous at those points.
         worst_gamma = max(worst_gamma, &
                           abs(f_gamma - vsigma(i))/max(abs(vsigma(i)), tiny(1.0_dp)))
      end do

      call check(error, n_live > 0, what//": no points above the density threshold")
      if (allocated(error)) return
      ! Non-vacuous: the second derivative's own sums must be alive too, or a
      ! future all-zero regression would sail through the identity above.
      call check(error, maxval(abs(a)) > 0.0_dp .and. maxval(abs(cc)) > 0.0_dp &
                 .and. maxval(abs(d2odr2)) > 0.0_dp .and. maxval(abs(d2kdr2)) > 0.0_dp, &
                 what//": a hessian intermediate is identically zero")
      if (allocated(error)) return
      if (worst_rho > thr .or. worst_gamma > thr_rel) then
         call test_failed(error, what//": rebuilt potential disagrees with vrho/vsigma")
         print '(a,es24.16)', "  worst f_rho   deviation          ", worst_rho
         print '(a,es24.16)', "  worst f_gamma relative deviation ", worst_gamma
      end if
   end subroutine rebuild_check

   subroutine test_rebuild_synthetic(error)
      !! The rebuild identity on the five reference points
      type(error_type), allocatable, intent(out) :: error
      call rebuild_check(B_VV, C_VV, COORDS, RHO, SIGMA, WEIGHTS, 1.0e-16_dp, &
                         1.0e-13_dp, "synthetic", error)
   end subroutine test_rebuild_synthetic

   subroutine test_rebuild_real_grid(error)
      !! The same identity on b97m-v's own NLC grid over a converged density
      !!
      !! The five synthetic points cannot reach what a molecular grid reaches
      !! every time: densities across ten orders of magnitude, near-zero
      !! gradients at the nuclei, and points that fall below the threshold
      !! mid-grid. This is the grid and density the Hessian will actually run
      !! on, so the conventions are pinned where they will be used.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: WATER_Z(3) = [8, 1, 1]
      character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
      real(dp), parameter :: WATER(3, 3) = reshape([ &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 1.4300_dp, 1.1075_dp, &
                                                   0.0000_dp, -1.4300_dp, 1.1075_dp], [3, 3])
      integer, parameter :: BLK = 256

      type(czt_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: rho(:), sigma(:)
      real(dp), allocatable :: ao(:, :), ao_grad(:, :, :)
      real(dp), allocatable :: rho_blk(:), rho_grad_blk(:, :)
      integer :: npts, g0, g1, nb, ig

      if (.not. xc_available()) return

      call build_czt_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "building water failed")
      if (allocated(error)) return
      call xc_context_create(mol, "b97m-v", ctx, err, level=3)
      if (err%has_error()) then
         call mol%destroy()
         call check(error, .false., "creating the b97m-v context failed")
         return
      end if
      call run_czt_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (.not. err%has_error()) call ensure_nlc_grid(ctx, mol, err)
      if (err%has_error()) then
         call mol%destroy()
         call check(error, .false., "the b97m-v reference state failed")
         return
      end if

      ! rho and sigma on the NLC grid, exactly as the potential path forms them
      npts = ctx%nlc_grid%n_points
      allocate (rho(npts), sigma(npts))
      do g0 = 1, npts, BLK
         g1 = min(g0 + BLK - 1, npts)
         nb = g1 - g0 + 1
         if (allocated(rho_blk)) deallocate (rho_blk, rho_grad_blk)
         allocate (rho_blk(nb), rho_grad_blk(nb, 3))
         call eval_ao_block(mol, ctx%nlc_grid%coords(:, g0:g1), ao, err, grad=ao_grad)
         if (err%has_error()) then
            call mol%destroy()
            call check(error, .false., "evaluating the density on the grid failed")
            return
         end if
         call eval_rho(ao, scf%density, rho_blk, ao_grad=ao_grad, rho_grad=rho_grad_blk)
         do ig = 1, nb
            rho(g0 + ig - 1) = rho_blk(ig)
            sigma(g0 + ig - 1) = rho_grad_blk(ig, 1)**2 + rho_grad_blk(ig, 2)**2 &
                                 + rho_grad_blk(ig, 3)**2
         end do
      end do
      call mol%destroy()

      call rebuild_check(ctx%nlc_b, ctx%nlc_c, ctx%nlc_grid%coords, rho, sigma, &
                         ctx%nlc_grid%weights, 1.0e-14_dp, 1.0e-13_dp, "real grid", error)
   end subroutine test_rebuild_real_grid

   subroutine run(exc, vrho, vsigma)
      real(dp), intent(out) :: exc(NPT), vrho(NPT), vsigma(NPT)
      call vv10_nlc(B_VV, C_VV, COORDS, RHO, SIGMA, &
                    COORDS, RHO, SIGMA, WEIGHTS, exc, vrho, vsigma)
   end subroutine run

   subroutine test_exc(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: exc(NPT), vrho(NPT), vsigma(NPT)

      call run(exc, vrho, vsigma)
      if (maxval(abs(exc - EXC_REF)) > THR) then
         call test_failed(error, "vv10 exc does not match")
         print '(a,es24.16)', "  max deviation ", maxval(abs(exc - EXC_REF))
         print '(3es24.16)', exc
      end if
   end subroutine test_exc

   subroutine test_potential(error)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: exc(NPT), vrho(NPT), vsigma(NPT)

      call run(exc, vrho, vsigma)
      if (maxval(abs(vrho - VRHO_REF)) > THR) then
         call test_failed(error, "vv10 vrho does not match")
         print '(a,es24.16)', "  max deviation ", maxval(abs(vrho - VRHO_REF))
         return
      end if
      ! vsigma is four orders smaller than vrho, so an absolute threshold
      ! would pass on a value of zero. Checked relatively.
      if (maxval(abs(vsigma - VSIGMA_REF)/max(abs(VSIGMA_REF), tiny(1.0_dp))) > 1.0e-10_dp) then
         call test_failed(error, "vv10 vsigma does not match")
         print '(3es24.16)', vsigma
      end if
   end subroutine test_potential

   subroutine test_energy(error)
      !! The quantity the SCF actually consumes, formed the way it will be
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: exc(NPT), vrho(NPT), vsigma(NPT), e_nl

      call run(exc, vrho, vsigma)
      e_nl = sum(WEIGHTS*RHO*exc)
      call check(error, e_nl, E_NL_REF, thr=THR)
      if (allocated(error)) return
      ! and it must not be zero, or every comparison above is vacuous
      call check(error, abs(e_nl) > 1.0e-6_dp, "the reference energy is not zero")
   end subroutine test_energy

   subroutine test_zero_gradient(error)
      !! sigma = 0 happens at every nucleus and at every local extremum
      !!
      !! `dw0_dsigma` divides by sigma while its numerator is proportional to
      !! sigma squared, so the limit is zero and the expression is 0/0. A NaN
      !! here would propagate into the Fock matrix and converge to nothing.
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: exc(NPT), vrho(NPT), vsigma(NPT)
      real(dp) :: flat(NPT)

      flat = 0.0_dp
      call vv10_nlc(B_VV, C_VV, COORDS, RHO, flat, &
                    COORDS, RHO, flat, WEIGHTS, exc, vrho, vsigma)
      if (any(exc /= exc) .or. any(vrho /= vrho) .or. any(vsigma /= vsigma)) then
         call test_failed(error, "a vanishing gradient produced NaN")
         return
      end if
      call check(error, maxval(abs(vsigma)), 0.0_dp, thr=0.0_dp)
   end subroutine test_zero_gradient

   subroutine test_threshold(error)
      !! A grid is mostly vacuum, and vacuum must contribute exactly nothing
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: exc(NPT), vrho(NPT), vsigma(NPT)
      real(dp) :: exc2(NPT + 2), vrho2(NPT + 2), vsigma2(NPT + 2)
      real(dp) :: c2(3, NPT + 2), r2(NPT + 2), s2(NPT + 2), w2(NPT + 2)

      call run(exc, vrho, vsigma)

      ! the same system with two empty points bolted on
      c2(:, 1:NPT) = COORDS; r2(1:NPT) = RHO; s2(1:NPT) = SIGMA; w2(1:NPT) = WEIGHTS
      c2(:, NPT + 1) = [5.0_dp, 5.0_dp, 5.0_dp]
      c2(:, NPT + 2) = [-6.0_dp, 2.0_dp, -4.0_dp]
      r2(NPT + 1:) = 0.1_dp*VV10_RHO_THRESHOLD
      s2(NPT + 1:) = 1.0e-20_dp
      w2(NPT + 1:) = 1.0_dp

      call vv10_nlc(B_VV, C_VV, c2, r2, s2, c2, r2, s2, w2, exc2, vrho2, vsigma2)

      if (maxval(abs(exc2(1:NPT) - exc)) > 0.0_dp) then
         call test_failed(error, "below-threshold points changed the result")
         return
      end if
      call check(error, maxval(abs(exc2(NPT + 1:))), 0.0_dp, thr=0.0_dp)
   end subroutine test_threshold

end module test_mqc_vv10

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_vv10, only: collect_vv10
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_vv10", collect_vv10)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
