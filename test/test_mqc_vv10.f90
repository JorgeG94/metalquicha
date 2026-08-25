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
   use mqc_libcint_vv10, only: vv10_nlc, VV10_RHO_THRESHOLD
   implicit none
   private

   public :: collect_vv10

   integer, parameter :: NPT = 5

   !> wB97X-V's parameters, as libxc reports them through xc_f03_nlc_coef
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

   !> The reference carries fourteen digits, so this is close to what double
   !> precision can express over a sum of this length. A looser threshold would
   !> accept a dropped term small enough to matter later.
   real(dp), parameter :: THR = 1.0e-12_dp

contains

   subroutine collect_vv10(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("vv10 energy density matches the reference", test_exc), &
                  new_unittest("vv10 potential matches the reference", test_potential), &
                  new_unittest("vv10 integrates to the reference energy", test_energy), &
                  new_unittest("vv10 survives a vanishing gradient", test_zero_gradient), &
                  new_unittest("vv10 ignores points below threshold", test_threshold) &
                  ]
   end subroutine collect_vv10

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
      !> vsigma is four orders smaller than vrho, so an absolute threshold
      !> would pass on a value of zero. Checked relatively.
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
      !> and it must not be zero, or every comparison above is vacuous
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

      !> the same system with two empty points bolted on
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
