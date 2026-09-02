!! The exchange-correlation Hessian against differences of its own gradient
module test_mqc_xc_hessian
   !! One test, and it is the one that decides whether the algebra is right.
   !!
   !! `xc_hessian` is `d2 E_xc / dR dR` with the grid held fixed, and
   !! `xc_gradient_fixed_grid` is the first derivative under the same
   !! assumption. Differencing the second against the first is therefore an
   !! exact comparison rather than an approximate one: any disagreement beyond
   !! the step error is a mistake in the second-derivative expression, which is
   !! where the mistakes are.
   !!
   !! Differenced against the *physical* `xc_gradient` this would disagree by
   !! the grid-response term it deliberately omits -- around 1e-4 on a gradient,
   !! which is small enough to be mistaken for a loose tolerance and large
   !! enough to hide a sign error. Hence the separate first derivative.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available, &
                             xc_add_potential, vv10_add_potential, &
                             vv10_kernel_apply, xc_grid_kernel_quantities, &
                             KERNEL_RHO_FLOOR, xc_kernel_apply, xc_kernel2_apply
   use mqc_libcint_hessian, only: ks_hessian
   use mqc_libcint_gradient, only: libcint_scf_gradient, vv10_gradient_fixed_grid, &
                                   xc_potential_gradient
   use mqc_libcint_xc_hessian, only: xc_hessian, xc_gradient_fixed_grid, &
                                     xc_potential_deriv, vv10_hessian, &
                                     vv10_potential_deriv, xc_potential_hessian, &
                                     xc_kernel_deriv
   implicit none
   private

   public :: collect_mqc_xc_hessian

   ! Water, bent, in Bohr. Small enough that a nine-by-nine finite-difference
   ! comparison is a second and large enough to have every atom pair in it.
   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape([ &
                                                0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                0.0000_dp, 1.4300_dp, 1.1075_dp, &
                                                0.0000_dp, -1.4300_dp, 1.1075_dp], [3, 3])

   ! The water dimer `test_vxc_potential_hessian_dimer` introduced, at module
   ! scope for the kernel tests that need the same things it was chosen for:
   ! d functions on both oxygens and no symmetry plane. Planar water at STO-3G
   ! zeroes every entry mixing an in-plane component with an out-of-plane one
   ! by reflection, which is exactly where a wrong term hides.
   integer, parameter :: DIMER_Z6(6) = [8, 1, 1, 8, 1, 1]
   character(len=2), parameter :: DIMER_SYM6(6) = ["O ", "H ", "H ", "O ", "H ", "H "]
   real(dp), parameter :: DIMER_XYZ6(3, 6) = reshape([ &
                                                     0.8974_dp, -1.285111_dp, 1.375674_dp, &
                                                     0.93366_dp, -1.620249_dp, 0.461291_dp, &
                                                     1.538424_dp, -0.56327_dp, 1.325981_dp, &
                                                     -1.559218_dp, -0.241778_dp, 1.423474_dp, &
                                                     -2.064102_dp, -0.501966_dp, 2.197464_dp, &
                                                     -0.677656_dp, -0.678646_dp, 1.525237_dp], &
                                                     [3, 6])

contains

   subroutine collect_mqc_xc_hessian(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("xc_hessian_differences_its_own_gradient", test_against_fd), &
                  new_unittest("xc_hessian_mgga_differences_its_gradient", test_gga_against_fd), &
                  new_unittest("xc_hessian_on_d_functions_differences_its_gradient", &
                               test_d_functions), &
                  new_unittest("vv10_fixed_grid_gradient_differences_the_energy", &
                               test_vv10_fixed_grid), &
                  new_unittest("vv10_hessian_differences_the_fixed_grid_gradient", &
                               test_vv10_hessian), &
                  new_unittest("xc_potential_derivative_matches_differences", test_vxc_deriv), &
                  new_unittest("vv10_potential_derivative_matches_differences", &
                               test_vv10_vxc_deriv), &
                  new_unittest("vv10_kernel_differences_the_potential_in_density", &
                               test_vv10_kernel_apply), &
                  new_unittest("third_derivative_differences_the_kernel_in_density", &
                               test_kxc_against_fxc), &
                  new_unittest("fixed_grid_potential_gradient_holds_its_grid", &
                               test_fixed_grid_potential_gradient), &
                  new_unittest("potential_hessian_on_a_dimer_with_d_functions", &
                               test_vxc_potential_hessian_dimer), &
                  new_unittest("potential_hessian_differences_the_potential_gradient", &
                               test_vxc_potential_hessian), &
                  new_unittest("kernel_derivative_differences_the_kernel_apply", &
                               test_kernel_deriv), &
                  new_unittest("kernel_derivative_on_a_dimer_with_d_functions", &
                               test_kernel_deriv_dimer), &
                  new_unittest("second_kernel_differences_the_kernel_in_density", &
                               test_kernel2), &
                  new_unittest("second_kernel_on_a_dimer_with_d_functions", &
                               test_kernel2_dimer), &
                  new_unittest("ks_hessian_differences_the_dft_gradient", test_ks_end_to_end) &
                  ]
   end subroutine collect_mqc_xc_hessian

   subroutine xc_at(ctx, coords, gradient, hess, density, err, ok, basis)
      !! A derivative at `coords`, on a grid that was built somewhere else
      !!
      !! Two things are held fixed and both matter. **The density matrix**,
      !! because `xc_hessian` is the explicit term -- what the energy does when
      !! the nuclei move and the density does not -- so re-converging at each
      !! displacement would add the whole coupled-perturbed contribution to one
      !! side only.
      !!
      !! **The grid**, because `xc_context_create` builds an atom-centred one
      !! and would therefore move it with the nuclei. The Hessian omits the grid
      !! response deliberately; a finite difference over rebuilt grids does not,
      !! and the two then disagree by that term. It is about 1.6e-4 here, which
      !! is both too small to look like an error and too large to be step noise
      !! -- precisely the size this repository already documents for the same
      !! term in a gradient. Passing one context in and moving only `mol` is
      !! what makes the comparison exact instead of approximate.
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: coords(3, 3)
      character(len=*), intent(in), optional :: basis
      real(dp), intent(out), optional :: gradient(3, 3)
      real(dp), intent(out), optional :: hess(3, 3, 3, 3)
      real(dp), intent(in) :: density(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      character(len=32) :: use_basis

      ok = .false.
      use_basis = "sto-3g"
      if (present(basis)) use_basis = basis
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, trim(use_basis), mol, err)
      if (err%has_error()) return
      if (present(gradient)) then
         gradient = 0.0_dp
         call xc_gradient_fixed_grid(ctx, mol, density, gradient, err)
      end if
      if (present(hess)) then
         hess = 0.0_dp
         call xc_hessian(ctx, mol, density, hess, err)
      end if
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine xc_at

   subroutine reference_state(functional, ctx, density, err, ok, basis)
      !! One converged reference: the grid and the density everything below uses
      !!
      !! `allow_half` is set here and in the two helpers below because several
      !! callers name a libxc exchange half on purpose -- `lda_x`, `gga_x_pbe`,
      !! `mgga_x_tpss` -- to get one rung's kernel with no correlation term beside
      !! it, so that a broken second derivative of exchange cannot be masked by a
      !! working one of correlation. A deck naming a half is refused, which is
      !! what that flag defaults to; asking for one here is deliberate, and the
      !! argument is how it says so.
      character(len=*), intent(in) :: functional
      character(len=*), intent(in), optional :: basis
         !! STO-3G unless asked otherwise, which is what every caller here
         !! wanted until d functions needed covering.
      type(xc_context_t), intent(out) :: ctx
      real(dp), allocatable, intent(out) :: density(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      character(len=32) :: use_basis

      ok = .false.
      use_basis = "sto-3g"
      if (present(basis)) use_basis = basis
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, trim(use_basis), mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, functional, ctx, err, level=3, allow_half=.true.)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (.not. err%has_error()) density = scf%density
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine reference_state

   subroutine test_against_fd(error)
      !! Every one of the nine by nine entries, against a central difference
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 1.0e-5_dp
         !! Ten times the step error of a central difference at this `H`, which
         !! measures about 1e-6 here, and far below anything a real mistake
         !! produces: leaving the density free to relax was 0.19 on this entry,
         !! and letting the grid move with the nuclei was 1.6e-4. The band
         !! between those and step noise is wide, so the tolerance is not a
         !! judgement call.
      real(dp) :: hess(3, 3, 3, 3), plus(3, 3), minus(3, 3), shifted(3, 3)
      real(dp) :: fd
      real(dp), allocatable :: dens(:, :)
      type(xc_context_t) :: ctx
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b

      if (.not. xc_available()) return

      call reference_state("lda_x", ctx, dens, err, ok)
      call check(error, ok, "the reference Kohn-Sham state failed")
      if (allocated(error)) return

      call xc_at(ctx, WATER, hess=hess, density=dens, err=err, ok=ok)
      call check(error, ok, "the analytic exchange-correlation Hessian failed")
      if (allocated(error)) return

      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call xc_at(ctx, shifted, gradient=plus, density=dens, err=err, ok=ok)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call xc_at(ctx, shifted, gradient=minus, density=dens, err=err, ok=ok)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            do ja = 1, 3
               do b = 1, 3
                  fd = (plus(b, ja) - minus(b, ja))/(2.0_dp*H)
                  call check(error, hess(a, b, ia, ja), fd, thr=TOL, &
                             more="exchange-correlation Hessian entry disagrees with "// &
                             "a difference of its own gradient")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_against_fd

   subroutine test_d_functions(error)
      !! The same exact comparison on a basis that has d functions
      !!
      !! Every check in this file ran on STO-3G, which is s and p only, so the
      !! angular part of a second derivative was never exercised past l = 1 and
      !! a term right for s and p and wrong for d would have passed all of it.
      !!
      !! **A meta-GGA is deliberately not here**, and the reason looks exactly
      !! like a bug if you meet it cold. On a level-3 grid `tpss`/cc-pVDZ
      !! appears to disagree with a difference of its own gradient by 4e-2 on
      !! the oxygen diagonal block. It does not. Swept over five step sizes
      !! from 1e-2 down to 6.25e-4 the difference returns 0.589, 0.792, 0.740,
      !! 0.768 and 0.748 -- scattering by 0.03 and converging on nothing, where
      !! a central difference of a smooth function converges as h^2. A second
      !! difference of the exchange-correlation *energy* scatters the same way,
      !! so it is the underlying function that is not smooth at this grid and
      !! not either derivative. M06-L is worse: 4.04, -0.30, 0.67, 0.53, 0.63.
      !!
      !! The tau channel is the most grid-sensitive of the three and d
      !! functions make it worse, so level 3 is simply too coarse to
      !! differentiate twice -- the Hessian itself is not grid-converged until
      !! level 7, where level 3 is 5.7e-3 out. At levels 7 and 9 the sweep does
      !! converge, and converges to this Hessian. Running it here at level 7
      !! would cost eighteen fine-grid gradients for a check the validation
      !! suite already makes better, against `pyscf.hessian.rks` element by
      !! element.
      type(error_type), allocatable, intent(out) :: error

      if (.not. xc_available()) return

      call d_function_fd_for("lda_x", error)
      if (allocated(error)) return
      call d_function_fd_for("pbe", error)
   end subroutine test_d_functions

   subroutine d_function_fd_for(functional, error)
      !! One functional on cc-pVDZ, `xc_hessian` against a difference of
      !! `xc_gradient_fixed_grid`
      character(len=*), intent(in) :: functional
      type(error_type), allocatable, intent(out) :: error

      character(len=*), parameter :: BASIS = "cc-pvdz"
      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 5.0e-4_dp
      real(dp) :: hess(3, 3, 3, 3), plus(3, 3), minus(3, 3), shifted(3, 3)
      real(dp) :: fd
      real(dp), allocatable :: dens(:, :)
      type(xc_context_t) :: ctx
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b

      call reference_state(functional, ctx, dens, err, ok, basis=BASIS)
      call check(error, ok, "the reference Kohn-Sham state failed for "//functional)
      if (allocated(error)) return

      call xc_at(ctx, WATER, hess=hess, density=dens, err=err, ok=ok, basis=BASIS)
      call check(error, ok, "the analytic Hessian failed for "//functional)
      if (allocated(error)) return

      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call xc_at(ctx, shifted, gradient=plus, density=dens, err=err, ok=ok, basis=BASIS)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call xc_at(ctx, shifted, gradient=minus, density=dens, err=err, ok=ok, basis=BASIS)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            do ja = 1, 3
               do b = 1, 3
                  fd = (plus(b, ja) - minus(b, ja))/(2.0_dp*H)
                  call check(error, hess(a, b, ia, ja), fd, thr=TOL, &
                             more="exchange-correlation Hessian entry disagrees with "// &
                             "a difference of its own gradient for "//functional//"/cc-pVDZ")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine d_function_fd_for

   subroutine test_gga_against_fd(error)
      !! The same comparison for a GGA, where the sigma channel is
      !!
      !! Worth its own test rather than a parameter on the first: a GGA adds
      !! five terms the LDA does not have, and the one that consumes third
      !! derivatives of the basis functions appears only on the diagonal atom
      !! block. An error confined there would leave every off-diagonal entry
      !! correct, so a test that stopped at the first disagreement on an LDA
      !! would never reach it.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 5.0e-4_dp
      real(dp) :: hess(3, 3, 3, 3), plus(3, 3), minus(3, 3), shifted(3, 3)
      real(dp) :: fd
      real(dp), allocatable :: dens(:, :)
      type(xc_context_t) :: ctx
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b

      if (.not. xc_available()) return

      call reference_state("mgga_x_tpss", ctx, dens, err, ok)
      call check(error, ok, "the reference Kohn-Sham state failed")
      if (allocated(error)) return

      call xc_at(ctx, WATER, hess=hess, density=dens, err=err, ok=ok)
      call check(error, ok, "the analytic meta-GGA exchange-correlation Hessian failed")
      if (allocated(error)) return

      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call xc_at(ctx, shifted, gradient=plus, density=dens, err=err, ok=ok)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call xc_at(ctx, shifted, gradient=minus, density=dens, err=err, ok=ok)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            do ja = 1, 3
               do b = 1, 3
                  fd = (plus(b, ja) - minus(b, ja))/(2.0_dp*H)
                  call check(error, hess(a, b, ia, ja), fd, thr=TOL, &
                             more="meta-GGA exchange-correlation Hessian entry "// &
                             "disagrees with a difference of its own gradient")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_gga_against_fd

   subroutine test_vv10_fixed_grid(error)
      !! The fixed-grid VV10 gradient against a difference of the fixed-grid energy
      !!
      !! This is scaffolding, and it is the scaffolding a VV10 Hessian cannot
      !! be built without. `vv10_gradient_fixed_grid` is what that Hessian will
      !! be differenced against, so an error in it would be inherited by every
      !! check made with it and would read as agreement.
      !!
      !! Both grids and the density are held fixed, by passing one context in
      !! and moving only the molecule -- the same trick `xc_at` uses, for the
      !! same reason. The semilocal term rides along rather than being isolated:
      !! it is validated already, so a failure here localises to the non-local
      !! part anyway, and `xc_add_potential` reports one energy for both.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp) :: shifted(3, 3), fd, worst, ep, em
      real(dp) :: gradient(3, 3)
      real(dp), allocatable :: density(:, :)
      type(xc_context_t) :: ctx
      type(error_t) :: err
      logical :: ok
      integer :: ia, a

      if (.not. xc_available()) return

      call reference_state("b97m-v", ctx, density, err, ok)
      call check(error, ok, "the b97m-v reference failed")
      if (allocated(error)) return

      gradient = 0.0_dp
      call xc_at(ctx, WATER, gradient=gradient, density=density, err=err, ok=ok)
      call check(error, ok, "the fixed-grid semilocal gradient failed")
      if (allocated(error)) return
      call vv10_at(ctx, WATER, gradient, density, err, ok)
      call check(error, ok, "the fixed-grid VV10 gradient failed")
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call exc_at(ctx, shifted, density, ep, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced exchange-correlation energy failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call exc_at(ctx, shifted, density, em, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced exchange-correlation energy failed")
               return
            end if
            fd = (ep - em)/(2.0_dp*H)
            worst = max(worst, abs(gradient(a, ia) - fd))
         end do
      end do

      ! Exact rather than approximate: nothing is omitted on either side, so
      ! what is left is the step error. Differencing the *physical* VV10
      ! gradient here instead would miss by 6e-3, the size of the terms this
      ! deliberately drops.
      call check(error, worst < 1.0e-6_dp, &
                 "the fixed-grid VV10 gradient disagrees with a difference of the energy")
   end subroutine test_vv10_fixed_grid

   subroutine test_vv10_hessian(error)
      !! The VV10 Hessian against differences of the fixed-grid VV10 gradient
      !!
      !! The comparison is exact in the sense that matters: both sides hold the
      !! density matrix and both grids fixed by passing one context in, and
      !! both omit the grid response, so nothing is approximated on one side
      !! only and what remains is the step error of the central difference.
      !! The VV10 term is isolated -- `vv10_hessian` against differences of
      !! `vv10_gradient_fixed_grid` alone -- because the semilocal term has its
      !! own test above and a failure here should localise to the non-local
      !! algebra, not to whichever term happens to be larger.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 1.0e-7_dp
         !! Six times the measured worst disagreement, 1.63e-8, which is the
         !! step error and nothing else: it scales exactly as H^2 over an
         !! eightfold range of steps (2.6e-7 at 1e-3, 1.04e-6 at 2e-3), which
         !! a missing term cannot do. The nearest real mistake is six orders
         !! away -- dropping the kernel's diagonal, the term where omega and
         !! kappa feel their own point's density, leaves every off-diagonal
         !! pair term intact and still misses by 4.6e-2 here.
      real(dp) :: hess(3, 3, 3, 3), plus(3, 3), minus(3, 3), shifted(3, 3)
      real(dp) :: fd, worst
      real(dp), allocatable :: density(:, :)
      type(xc_context_t) :: ctx
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b

      if (.not. xc_available()) return

      call reference_state("b97m-v", ctx, density, err, ok)
      call check(error, ok, "the b97m-v reference failed")
      if (allocated(error)) return

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      hess = 0.0_dp
      call vv10_hessian(ctx, mol, density, hess, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the analytic VV10 Hessian failed")
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            plus = 0.0_dp
            call vv10_at(ctx, shifted, plus, density, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced fixed-grid VV10 gradient failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            minus = 0.0_dp
            call vv10_at(ctx, shifted, minus, density, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced fixed-grid VV10 gradient failed")
               return
            end if
            do ja = 1, 3
               do b = 1, 3
                  fd = (plus(b, ja) - minus(b, ja))/(2.0_dp*H)
                  worst = max(worst, abs(hess(a, b, ia, ja) - fd))
               end do
            end do
         end do
      end do

      call check(error, worst < TOL, &
                 "the VV10 Hessian disagrees with a difference of the "// &
                 "fixed-grid VV10 gradient")
   end subroutine test_vv10_hessian

   subroutine vv10_at(ctx, coords, gradient, density, err, ok)
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: coords(3, 3)
      real(dp), intent(inout) :: gradient(3, 3)
      real(dp), intent(in) :: density(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call vv10_gradient_fixed_grid(ctx, mol, density, gradient, err)
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine vv10_at

   subroutine exc_at(ctx, coords, density, e_xc, err, ok)
      !! E_xc, non-local term included, on the grid the context already holds
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: coords(3, 3)
      real(dp), intent(in) :: density(:, :)
      real(dp), intent(out) :: e_xc
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: v(:, :)
      real(dp) :: n_elec

      ok = .false.
      e_xc = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      allocate (v(mol%nao, mol%nao))
      v = 0.0_dp
      call xc_add_potential(ctx, mol, density, v, e_xc, n_elec, err)
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine exc_at

   subroutine test_vxc_deriv(error)
      !! The exchange-correlation potential's nuclear derivative, against
      !! central differences of the potential matrix itself
      !!
      !! This is the term the response half of a Kohn-Sham Hessian needs and
      !! that `h1_contract` cannot supply, because it knows only about
      !! integrals. Differenced against `xc_add_potential` -- the same matrix
      !! the SCF builds -- with the density matrix and the grid both held fixed,
      !! for the reasons the Hessian test above sets out.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp) :: shifted(3, 3), worst
      real(dp), allocatable :: dens(:, :), h1(:, :, :, :), vp(:, :), vm(:, :)
      type(xc_context_t) :: ctx
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, nao, u, v

      if (.not. xc_available()) return

      call vxc_deriv_for("lda_x", error)
      if (allocated(error)) return
      call vxc_deriv_for("gga_x_pbe", error)
      if (allocated(error)) return
      call vxc_deriv_for("pbe", error)
      if (allocated(error)) return
      call vxc_deriv_for("tpss", error)
   end subroutine test_vxc_deriv

   subroutine vxc_deriv_for(functional, error)
      !! One functional, every element of dV/dR against a central difference
      character(len=*), intent(in) :: functional
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp) :: shifted(3, 3), worst
      real(dp), allocatable :: dens(:, :), h1(:, :, :, :), vp(:, :), vm(:, :)
      type(xc_context_t) :: ctx
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, nao, u, v

      call reference_state(functional, ctx, dens, err, ok)
      call check(error, ok, "the reference Kohn-Sham state failed")
      if (allocated(error)) return

      nao = size(dens, 1)
      allocate (h1(nao, nao, 3, size(WATER_Z)))
      h1 = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call xc_potential_deriv(ctx, mol, dens, h1, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the potential derivative failed")
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call vxc_at(ctx, shifted, dens, vp, err)
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call vxc_at(ctx, shifted, dens, vm, err)
            do v = 1, nao
               do u = 1, nao
                  worst = max(worst, abs(h1(u, v, a, ia) - (vp(u, v) - vm(u, v))/(2.0_dp*H)))
               end do
            end do
         end do
      end do

      call check(error, worst < 1.0e-5_dp, &
                 "the exchange-correlation potential derivative disagrees with "// &
                 "differences of the potential for "//functional)
   end subroutine vxc_deriv_for

   subroutine vxc_at(ctx, coords, dens, v, err)
      !! The exchange-correlation potential matrix at a displaced geometry, on
      !! the reference grid and against the reference density
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: coords(3, 3), dens(:, :)
      real(dp), allocatable, intent(out) :: v(:, :)
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      real(dp) :: e, n

      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (allocated(v)) deallocate (v)
      allocate (v(size(dens, 1), size(dens, 2)))
      v = 0.0_dp
      call xc_add_potential(ctx, mol, dens, v, e, n, err)
      call mol%destroy()
   end subroutine vxc_at

   subroutine test_vv10_vxc_deriv(error)
      !! The VV10 potential's nuclear derivative, against central differences
      !! of the VV10 potential matrix itself
      !!
      !! `test_vxc_deriv`'s comparison for the non-local term: the density
      !! matrix and both grids held fixed by passing one context in and moving
      !! only the molecule, so both sides omit the grid response and what
      !! remains is the step error. The VV10 potential is isolated by calling
      !! `vv10_add_potential` directly into a zeroed matrix -- the same routine
      !! that adds VV10's term to `xc_add_potential`'s output -- so the
      !! semilocal term, which has its own test above, never enters and a
      !! failure localises to the non-local algebra.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 6.0e-8_dp
         !! Six times the measured worst disagreement, 1.01e-8 at this step,
         !! which is the step error and nothing else: it scales exactly as H^2
         !! over a fourfold range of steps (2.52e-9 at 5e-4, 4.03e-8 at 2e-3,
         !! ratios 4.00 and 4.00), which a missing term cannot do.
      real(dp) :: shifted(3, 3), worst
      real(dp), allocatable :: dens(:, :), h1(:, :, :, :), vp(:, :), vm(:, :)
      type(xc_context_t) :: ctx
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, nao, u, v

      if (.not. xc_available()) return

      call reference_state("b97m-v", ctx, dens, err, ok)
      call check(error, ok, "the b97m-v reference failed")
      if (allocated(error)) return

      nao = size(dens, 1)
      allocate (h1(nao, nao, 3, size(WATER_Z)))
      h1 = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call vv10_potential_deriv(ctx, mol, dens, h1, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the VV10 potential derivative failed")
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call vv10_vxc_at(ctx, shifted, dens, vp, err)
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call vv10_vxc_at(ctx, shifted, dens, vm, err)
            do v = 1, nao
               do u = 1, nao
                  worst = max(worst, abs(h1(u, v, a, ia) - (vp(u, v) - vm(u, v))/(2.0_dp*H)))
               end do
            end do
         end do
      end do

      call check(error, worst < TOL, &
                 "the VV10 potential derivative disagrees with differences "// &
                 "of the VV10 potential")
   end subroutine test_vv10_vxc_deriv

   subroutine vv10_vxc_at(ctx, coords, dens, v, err)
      !! The VV10 potential matrix alone, at a displaced geometry, on the NLC
      !! grid the context already holds and against the reference density
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: coords(3, 3), dens(:, :)
      real(dp), allocatable, intent(out) :: v(:, :)
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      real(dp) :: e

      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (allocated(v)) deallocate (v)
      allocate (v(size(dens, 1), size(dens, 2)))
      v = 0.0_dp
      call vv10_add_potential(ctx, mol, dens, v, e, err)
      call mol%destroy()
   end subroutine vv10_vxc_at

   subroutine test_vv10_kernel_apply(error)
      !! The VV10 kernel applied to trial densities, against central
      !! differences of the VV10 potential matrix in the *density*
      !!
      !! The Phase 4 comparison: `vv10_kernel_apply` is `dV/dD` contracted
      !! against a trial `P`, so the check is `(V(D + h P) - V(D - h P))/2h`
      !! with the geometry, the context and hence both grids never moving at
      !! all -- only the density matrix handed to `vv10_add_potential` changes.
      !! The VV10 potential is isolated exactly as `test_vv10_vxc_deriv`
      !! isolates it, so the semilocal kernel, which has `xc_kernel_apply`'s
      !! own tests, never enters.
      !!
      !! Three trials in one batched call, deliberately unlike each other: the
      !! density itself, a smooth off-diagonal band, and a diagonally dominant
      !! matrix -- a kernel wrong in only one channel has to disagree with at
      !! least one of them.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 6.0e-10_dp
         !! Six times the measured worst disagreement, 8.88e-11 at this step,
         !! which is the step error and nothing else: it scales exactly as H^2
         !! over a fourfold range of steps (2.22e-11 at 5e-4, 3.55e-10 at 2e-3,
         !! ratios 4.00 and 4.00), which a missing term cannot do. A break test
         !! dropping the `f_gamma grad drho` term misses by 4.67e-4,
         !! step-independent -- six orders above this band.
      integer, parameter :: N_TRIAL = 3
      real(dp) :: worst
      real(dp), allocatable :: dens(:, :), trial(:, :, :), vk(:, :, :)
      real(dp), allocatable :: vp(:, :), vm(:, :)
      type(xc_context_t) :: ctx
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      logical :: ok
      integer :: nao, u, v, it

      if (.not. xc_available()) return

      call reference_state("b97m-v", ctx, dens, err, ok)
      call check(error, ok, "the b97m-v reference failed")
      if (allocated(error)) return

      nao = size(dens, 1)
      allocate (trial(nao, nao, N_TRIAL), vk(nao, nao, N_TRIAL))
      trial(:, :, 1) = dens
      do v = 1, nao
         do u = 1, nao
            trial(u, v, 2) = 0.1_dp*cos(real(u - v, dp))
            trial(u, v, 3) = 0.05_dp*cos(real(u + v, dp))
         end do
         trial(v, v, 3) = trial(v, v, 3) + 0.2_dp
      end do

      ! Zeroed here because `vv10_kernel_apply` accumulates, as
      ! `xc_kernel_apply` does -- the convention its callers already carry.
      vk = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call vv10_kernel_apply(ctx, mol, dens, trial, vk, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the VV10 kernel apply failed")
      if (allocated(error)) return

      worst = 0.0_dp
      do it = 1, N_TRIAL
         call vv10_vxc_at(ctx, WATER, dens + H*trial(:, :, it), vp, err)
         call vv10_vxc_at(ctx, WATER, dens - H*trial(:, :, it), vm, err)
         do v = 1, nao
            do u = 1, nao
               worst = max(worst, abs(vk(u, v, it) - (vp(u, v) - vm(u, v))/(2.0_dp*H)))
            end do
         end do
      end do
      call check(error,.not. err%has_error(), "the displaced-density potentials failed")
      if (allocated(error)) return

      call check(error, worst < TOL, &
                 "the VV10 kernel disagrees with differences of the VV10 "// &
                 "potential in the density")
   end subroutine test_vv10_kernel_apply

   subroutine test_kxc_against_fxc(error)
      !! The third functional derivative against central differences of the
      !! second, in the density -- the rung above `xc_kernel_apply`'s
      !!
      !! The double-hybrid Hessian is the first thing here to differentiate the
      !! *kernel*, so `g_xc` is new and has nothing above it to be checked
      !! against. It has something below it: `f_xc`, which the ladder in this
      !! file already trusts. Perturbing the density matrix by a trial `P` moves
      !! both of the functional's arguments,
      !!
      !!     d rho   = rho_P
      !!     d sigma = 2 grad rho . grad rho_P
      !!
      !! so the chain rule ties every one of libxc's four third-derivative
      !! channels to a difference of the second ones:
      !!
      !!     d(v2rho2)     = v3rho3        d rho + v3rho2sigma   d sigma
      !!     d(v2rhosigma) = v3rho2sigma   d rho + v3rhosigma2   d sigma
      !!     d(v2sigma2)   = v3rhosigma2   d rho + v3sigma3      d sigma
      !!
      !! All four appear, and each in two different equations, so a channel
      !! transposed with its neighbour -- the mistake this bookkeeping invites --
      !! cannot satisfy both. Both families are run: B88 exercises all four,
      !! SVWN only `v3rho3`, which is the one an LDA-only implementation would
      !! get right by accident while the sigma channels stayed zero.
      !!
      !! The grid never moves. This is a derivative in the density at fixed
      !! geometry, so nothing here depends on the quadrature responding.
      type(error_type), allocatable, intent(out) :: error

      call kxc_for("gga_x_b88", error)
      if (allocated(error)) return
      call kxc_for("lda_x", error)
   end subroutine test_kxc_against_fxc

   subroutine kxc_for(functional, error)
      !! One functional, all channels the family carries
      character(len=*), intent(in) :: functional
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 5.0e-8_dp
         !! Seven times the measured worst disagreement, 7.33e-09 for B88 and
         !! 1.45e-09 for the LDA at this step. That is the step error and
         !! nothing else: it scales exactly as H^2 -- 2.93e-08 at 2e-4 against
         !! 7.33e-09 at 1e-4, a ratio of 4.00 -- which a missing or misplaced
         !! term cannot do.
         !!
         !! A break test transposing `v3rho2sigma` with `v3rhosigma2`, which is
         !! the mistake this bookkeeping invites, misses by 2.0 -- eight orders
         !! above this band and step-independent.
         !!
         !! Relative rather than absolute, because these are third derivatives
         !! of a functional that diverges at the tail of every atomic grid.
      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), trial(:, :)
      real(dp), allocatable :: r0(:), g0(:, :), vr0(:), vs0(:), f0rr(:), f0rs(:), f0ss(:)
      real(dp), allocatable :: krrr(:), krrs(:), krss(:), ksss(:)
      real(dp), allocatable :: rp(:), gp(:, :), vrp(:), vsp(:), fprr(:), fprs(:), fpss(:)
      real(dp), allocatable :: rm(:), gm(:, :), vrm(:), vsm(:), fmrr(:), fmrs(:), fmss(:)
      real(dp) :: drho, dsig, want, got, scale, worst
      integer :: n, i, j, ig, npts

      ! Setup failures are failures of this test. Returning with `error`
      ! unallocated is what testdrive counts as a pass, so a molecule that would
      ! not build or an SCF that would not converge used to leave the case green
      ! having compared nothing.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule did not build: "// &
                 err%get_message())
      if (allocated(error)) return
      call xc_context_create(mol, functional, ctx, err, level=3, allow_half=.true.)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the functional did not resolve: "// &
                 err%get_message())
      if (allocated(error)) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the reference did not converge: "// &
                 err%get_message())
      if (allocated(error)) return
      dens = scf%density
      n = size(dens, 1)

      ! A trial that is not proportional to the density: a proportional one
      ! moves rho and sigma together and cannot separate the channels.
      allocate (trial(n, n))
      do j = 1, n
         do i = 1, n
            trial(i, j) = 0.05_dp/(1.0_dp + real(abs(i - j), dp)) + 0.01_dp*dens(i, j)
         end do
      end do
      trial = 0.5_dp*(trial + transpose(trial))

      call xc_grid_kernel_quantities(ctx, mol, dens, r0, g0, vr0, vs0, f0rr, f0rs, f0ss, &
                                     err, grrr=krrr, grrs=krrs, grss=krss, gsss=ksss)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the kernel evaluation failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call xc_grid_kernel_quantities(ctx, mol, dens + H*trial, rp, gp, vrp, vsp, &
                                     fprr, fprs, fpss, err)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the kernel evaluation failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call xc_grid_kernel_quantities(ctx, mol, dens - H*trial, rm, gm, vrm, vsm, &
                                     fmrr, fmrs, fmss, err)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the kernel evaluation failed: "// &
                 err%get_message())
      if (allocated(error)) return

      npts = size(r0)
      worst = 0.0_dp
      do ig = 1, npts
         ! Skip the floored tail. Both orders are zeroed below the kernel's
         ! density floor, so differencing across that boundary compares a
         ! discontinuity against a derivative -- a comparison with no content.
         ! Ten times the floor keeps the whole stencil clear of it.
         if (min(r0(ig), rp(ig), rm(ig)) < 10.0_dp*KERNEL_RHO_FLOOR) cycle
         drho = (rp(ig) - rm(ig))/(2.0_dp*H)
         dsig = 0.0_dp
         do i = 1, 3
            dsig = dsig + 2.0_dp*g0(ig, i)*(gp(ig, i) - gm(ig, i))/(2.0_dp*H)
         end do

         ! Relative to the size of the terms being compared: at the tail of an
         ! atomic grid these are enormous and their difference is meaningless.
         want = (fprr(ig) - fmrr(ig))/(2.0_dp*H)
         got = krrr(ig)*drho + krrs(ig)*dsig
         scale = max(1.0_dp, abs(want), abs(got))
         worst = max(worst, abs(want - got)/scale)

         want = (fprs(ig) - fmrs(ig))/(2.0_dp*H)
         got = krrs(ig)*drho + krss(ig)*dsig
         scale = max(1.0_dp, abs(want), abs(got))
         worst = max(worst, abs(want - got)/scale)

         want = (fpss(ig) - fmss(ig))/(2.0_dp*H)
         got = krss(ig)*drho + ksss(ig)*dsig
         scale = max(1.0_dp, abs(want), abs(got))
         worst = max(worst, abs(want - got)/scale)
      end do

      call mol%destroy()
      call check(error, worst < TOL, trim(functional)// &
                 ": the third derivative should difference the second")
   end subroutine kxc_for

   subroutine test_fixed_grid_potential_gradient(error)
      !! `xc_potential_gradient(fixed_grid=.true.)` against the fixed-grid
      !! derivative this file already pins
      !!
      !! `P` does not respond, so `d Tr(P V) = Tr(P dV)` exactly, and
      !! `xc_potential_deriv` is that `dV/dR` with the grid held fixed. The two
      !! are the same number by construction, which is what makes this a check
      !! on the flag rather than on the physics.
      !!
      !! It exists because the flag shipped broken. Holding the grid means
      !! dropping *two* things -- the partition weights responding, and the
      !! points themselves travelling with their owning atom -- and the first
      !! version dropped only the first: `accumulate_channel` and
      !! `accumulate_gga_channel` take a `moving` argument that the four call
      !! sites never forwarded, though `vv10_gradient_core` forwards it and was
      !! the model to follow. What let it through review is that the *physical*
      !! path was untouched and every validation case still passed. Nothing
      !! exercised the flag at all.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: TOL = 1.0e-10_dp
         !! Two evaluations of one quantity, so this is round-off. There is no
         !! step here to carry a step error.
      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), pmat(:, :), grad(:, :), h1(:, :, :, :)
      real(dp) :: want, worst
      integer :: n, ia, c, u, v

      ! Setup failures are failures of this test. Returning with `error`
      ! unallocated is what testdrive counts as a pass, so a molecule that would
      ! not build or an SCF that would not converge used to leave the case green
      ! having compared nothing.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule did not build: "// &
                 err%get_message())
      if (allocated(error)) return
      call xc_context_create(mol, "gga_x_b88", ctx, err, level=3, allow_half=.true.)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the functional did not resolve: "// &
                 err%get_message())
      if (allocated(error)) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the reference did not converge: "// &
                 err%get_message())
      if (allocated(error)) return
      dens = scf%density
      n = size(dens, 1)

      ! Not the reference density. A second density that is not proportional to
      ! the first is what separates the two terms being dropped.
      allocate (pmat(n, n))
      do v = 1, n
         do u = 1, n
            pmat(u, v) = 0.03_dp/(1.0_dp + real(abs(u - v), dp)) + 0.01_dp*dens(u, v)
         end do
      end do
      pmat = 0.5_dp*(pmat + transpose(pmat))

      allocate (grad(3, 3), h1(n, n, 3, 3))
      grad = 0.0_dp
      h1 = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "the molecule did not build: "// &
                 err%get_message())
      if (allocated(error)) return
      call xc_potential_gradient(ctx, mol, dens, pmat, grad, err, fixed_grid=.true.)
      if (.not. err%has_error()) call xc_potential_deriv(ctx, mol, dens, h1, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do c = 1, 3
            want = 0.0_dp
            do v = 1, n
               do u = 1, n
                  want = want + pmat(u, v)*h1(u, v, c, ia)
               end do
            end do
            worst = max(worst, abs(grad(c, ia) - want))
         end do
      end do

      call check(error, worst < TOL, &
                 "the fixed-grid potential gradient should equal Tr(P dV/dR)")
   end subroutine test_fixed_grid_potential_gradient

   subroutine test_vxc_potential_hessian(error)
      !! The linear form's second derivative against central differences of
      !! its own fixed-grid first derivative
      !!
      !! `xc_potential_hessian` is `d2/dRdR Tr(P V_xc[D])` with the grid and
      !! both matrices fixed, and `Tr(P dV/dX)` -- `xc_potential_deriv`
      !! contracted with the fixed `P` -- is its first derivative under
      !! exactly the same approximation, so the comparison is exact up to step
      !! error: the construction of `xc_hessian` against
      !! `xc_gradient_fixed_grid`, one linear form over, differencing a rung
      !! this file already pins. Why the reference is not
      !! `xc_potential_gradient` with `fixed_grid` is documented on the helper
      !! below, with the numbers. Both families run, and B88 is the one that
      !! decides: an LDA gets the rho channel right while every sigma channel
      !! stays zero, so an implementation wrong in nothing but sigma would
      !! pass it.
      type(error_type), allocatable, intent(out) :: error

      if (.not. xc_available()) return

      call vxc_potential_hessian_for("lda_x", error)
      if (allocated(error)) return
      call vxc_potential_hessian_for("gga_x_b88", error)
   end subroutine test_vxc_potential_hessian

   subroutine vxc_potential_hessian_for(functional, error)
      !! One functional, all eighty-one entries against a central difference
      character(len=*), intent(in) :: functional
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 1.6e-8_dp
         !! Seven times the measured worst disagreement, 2.23e-9 for B88 at
         !! this step and 2.17e-10 for the LDA. The LDA residual is the step
         !! error and nothing else: 5.45e-11 at 5e-5, 8.67e-10 at 2e-4 and
         !! 3.47e-9 at 4e-4 -- ratios 3.98, 3.99, 4.00, exactly H^2. B88
         !! follows the same H^2 law under a ~2.2e-9 floor that moves
         !! *non-monotonically* with the step (2.27e-9 at 5e-5), which a wrong
         !! or missing term cannot do -- the break test dropping the g_xc
         !! group misses by 3.6521506, step-independent to nine digits over an
         !! eightfold range of steps. The floor is the quadrature tail: the
         !! displaced-geometry reference crosses KERNEL_RHO_FLOOR at points
         !! where B88's kernels are enormous, the discontinuity `kxc_for`
         !! excludes for the same reason, and nine orders below what a real
         !! mistake produces.
      real(dp) :: hess(3, 3, 3, 3), plus(3, 3), minus(3, 3), shifted(3, 3)
      real(dp) :: fd, worst
      real(dp), allocatable :: dens(:, :), trial(:, :)
      type(xc_context_t) :: ctx
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b, i, j, n

      call reference_state(functional, ctx, dens, err, ok)
      call check(error, ok, "the reference Kohn-Sham state failed for "//functional)
      if (allocated(error)) return

      ! The trial `kxc_for` uses, for the reason it gives: one proportional to
      ! the density moves rho and sigma together and cannot separate the
      ! channels. Symmetric and indefinite, like the relaxed density it stands
      ! in for.
      n = size(dens, 1)
      allocate (trial(n, n))
      do j = 1, n
         do i = 1, n
            trial(i, j) = 0.05_dp/(1.0_dp + real(abs(i - j), dp)) + 0.01_dp*dens(i, j)
         end do
      end do
      trial = 0.5_dp*(trial + transpose(trial))

      hess = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      call xc_potential_hessian(ctx, mol, dens, trial, hess, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), &
                 "the potential Hessian failed for "//functional)
      if (allocated(error)) return

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call vxc_potential_gradient_at(ctx, shifted, dens, trial, plus, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced potential gradient failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call vxc_potential_gradient_at(ctx, shifted, dens, trial, minus, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced potential gradient failed")
               return
            end if
            do ja = 1, 3
               do b = 1, 3
                  fd = (plus(b, ja) - minus(b, ja))/(2.0_dp*H)
                  worst = max(worst, abs(hess(a, b, ia, ja) - fd))
               end do
            end do
         end do
      end do

      call check(error, worst < TOL, &
                 "the potential Hessian disagrees with a difference of the "// &
                 "fixed-grid potential gradient for "//functional)
   end subroutine vxc_potential_hessian_for

   subroutine test_vxc_potential_hessian_dimer(error)
      !! The same rung on a water dimer in 6-31G(d): d functions, six atoms,
      !! and no symmetry to hide behind
      !!
      !! Water at STO-3G is planar and s,p only, and that combination masks
      !! things. Entries mixing an in-plane component with an out-of-plane one
      !! vanish by reflection there, which is exactly how the `fixed_grid` bug
      !! looked plausible for as long as it did -- it missed by 5.7e-02 on the
      !! in-plane components and agreed to step error on the rest. A dimer has
      !! no such plane, and 6-31G(d) puts Cartesian d functions on both oxygens,
      !! so the second derivatives of the basis reach angular momenta the
      !! monomer case never forms.
      !!
      !! Held to a looser band than the monomer for one reason: the quantity is
      !! larger, and this is an absolute comparison of entries whose own
      !! magnitudes are an order up.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: DZ(6) = [8, 1, 1, 8, 1, 1]
      character(len=2), parameter :: DSYM(6) = ["O ", "H ", "H ", "O ", "H ", "H "]
      real(dp), parameter :: DIMER(3, 6) = reshape([ &
                                                   0.8974_dp, -1.285111_dp, 1.375674_dp, &
                                                   0.93366_dp, -1.620249_dp, 0.461291_dp, &
                                                   1.538424_dp, -0.56327_dp, 1.325981_dp, &
                                                   -1.559218_dp, -0.241778_dp, 1.423474_dp, &
                                                   -2.064102_dp, -0.501966_dp, 2.197464_dp, &
                                                   -0.677656_dp, -0.678646_dp, 1.525237_dp], [3, 6])
      real(dp), parameter :: H = 1.0e-4_dp
      real(dp), parameter :: TOL = 2.0e-8_dp
         !! Six times the measured worst, 3.21e-09 across all 324 entries.
         !! Zeroing the `g_xc` group misses by 3.18 instead -- nine orders up,
         !! and on a system where no reflection plane can hide it.

      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), pmat(:, :), hess(:, :, :, :)
      real(dp), allocatable :: h1(:, :, :, :), coords(:, :)
      real(dp), allocatable :: gplus(:, :), gminus(:, :)
      real(dp) :: fd, worst
      integer :: n, ia, ja, a, b, u, v

      if (.not. xc_available()) return

      call build_libcint_molecule(DZ, DSYM, DIMER, "6-31g*", mol, err)
      call check(error,.not. err%has_error(), "the dimer did not build: "// &
                 err%get_message())
      if (allocated(error)) return
      call xc_context_create(mol, "gga_x_b88", ctx, err, level=3, allow_half=.true.)
      if (err%has_error()) call mol%destroy()
      call check(error,.not. err%has_error(), "the functional did not resolve: "// &
                 err%get_message())
      if (allocated(error)) return
      call run_libcint_rhf(mol, 20, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the dimer reference did not "// &
                 "converge: "//err%get_message())
      if (allocated(error)) return
      dens = scf%density
      n = size(dens, 1)

      allocate (pmat(n, n))
      do v = 1, n
         do u = 1, n
            pmat(u, v) = 0.02_dp/(1.0_dp + real(abs(u - v), dp)) + 0.01_dp*dens(u, v)
         end do
      end do
      pmat = 0.5_dp*(pmat + transpose(pmat))

      allocate (hess(3, 3, 6, 6))
      hess = 0.0_dp
      call build_libcint_molecule(DZ, DSYM, DIMER, "6-31g*", mol, err)
      call check(error,.not. err%has_error(), "the dimer did not rebuild: "// &
                 err%get_message())
      if (allocated(error)) return
      call xc_potential_hessian(ctx, mol, dens, pmat, hess, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      ! Every one of the 18 by 18 entries, not a sum rule over them. An
      ! aggregate would pass with errors that cancel across components, which is
      ! the failure mode a dimer was chosen to expose in the first place.
      allocate (coords(3, 6), h1(n, n, 3, 6), gplus(3, 6), gminus(3, 6))
      worst = 0.0_dp
      do ja = 1, 6
         do b = 1, 3
            coords = DIMER
            coords(b, ja) = DIMER(b, ja) + H
            call dimer_trace(DZ, DSYM, coords, ctx, dens, pmat, h1, gplus, err)
            if (err%has_error()) exit
            coords(b, ja) = DIMER(b, ja) - H
            call dimer_trace(DZ, DSYM, coords, ctx, dens, pmat, h1, gminus, err)
            if (err%has_error()) exit
            do ia = 1, 6
               do a = 1, 3
                  fd = (gplus(a, ia) - gminus(a, ia))/(2.0_dp*H)
                  worst = max(worst, abs(hess(a, b, ia, ja) - fd))
               end do
            end do
         end do
         if (err%has_error()) exit
      end do
      call check(error,.not. err%has_error(), err%get_message())
      if (allocated(error)) return

      call check(error, worst < TOL, &
                 "the potential Hessian should difference its gradient on a dimer")
   end subroutine test_vxc_potential_hessian_dimer

   subroutine dimer_trace(z, sym, coords, ctx, dens, pmat, h1, trace, err)
      !! Tr(P dV/dR) for every component, at one geometry
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: coords(:, :)
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: dens(:, :), pmat(:, :)
      real(dp), intent(inout) :: h1(:, :, :, :)
      real(dp), intent(out) :: trace(:, :)   !! (3, natm), Tr(P dV/dR) per component
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol
      integer :: ia, c, u, v
      real(dp) :: acc

      trace = 0.0_dp
      h1 = 0.0_dp
      call build_libcint_molecule(z, sym, coords, "6-31g*", mol, err)
      if (err%has_error()) return
      call xc_potential_deriv(ctx, mol, dens, h1, err)
      call mol%destroy()
      if (err%has_error()) return
      do ia = 1, size(z)
         do c = 1, 3
            acc = 0.0_dp
            do v = 1, size(dens, 1)
               do u = 1, size(dens, 1)
                  acc = acc + pmat(u, v)*h1(u, v, c, ia)
               end do
            end do
            trace(c, ia) = acc
         end do
      end do
   end subroutine dimer_trace

   subroutine vxc_potential_gradient_at(ctx, coords, dens, pmat, gradient, err, ok)
      !! The fixed-grid first derivative of `Tr(P V_xc[D])` at displaced
      !! coordinates: `xc_potential_deriv` contracted with `P`
      !!
      !! `P` is fixed, so `d/dX Tr(P V) = Tr(P dV/dX)` exactly, and
      !! `xc_potential_deriv` is the rung below this one on the ladder --
      !! `test_vxc_deriv` pins it against differences of the potential matrix
      !! itself. It carries no grid terms at all, which is the approximation
      !! the Hessian above makes.
      !!
      !! `xc_potential_gradient(fixed_grid=.true.)` is the same number and
      !! would serve, but this route shares less machinery with the routine
      !! under test, so a common-mode error has one fewer place to hide.
      !!
      !! That flag was broken when this test was written: its four
      !! `accumulate_channel` deposits did not forward `moving`, so it dropped
      !! the partition-weight response while *keeping* the grid-point-motion
      !! term -- against a scalar difference of `Tr(P V)` on the frozen grid it
      !! missed by up to 5.7e-2 here, in-plane components only, the
      !! out-of-plane ones agreeing to step error, which is what let it look
      !! plausible. It is fixed now, and `test_fixed_grid_potential_gradient`
      !! pins it. The history is kept because the shape of that bug is the
      !! reason this test reaches for the independent route.
      type(xc_context_t), intent(inout) :: ctx
      real(dp), intent(in) :: coords(3, 3), dens(:, :), pmat(:, :)
      real(dp), intent(out) :: gradient(3, 3)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      real(dp), allocatable :: h1(:, :, :, :)
      integer :: ia, a

      ok = .false.
      gradient = 0.0_dp
      allocate (h1(size(dens, 1), size(dens, 2), 3, size(WATER_Z)))
      h1 = 0.0_dp
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call xc_potential_deriv(ctx, mol, dens, h1, err)
      call mol%destroy()
      if (err%has_error()) return
      do ia = 1, size(WATER_Z)
         do a = 1, 3
            gradient(a, ia) = sum(pmat*h1(:, :, a, ia))
         end do
      end do
      ok = .true.
   end subroutine vxc_potential_gradient_at

   subroutine test_ks_end_to_end(error)
      !! The assembled Kohn-Sham Hessian against differences of the analytic
      !! density-functional gradient
      !!
      !! Every piece above is verified in isolation; this is the one that
      !! catches them being combined wrongly. Two independent signals are read:
      !! agreement with the finite difference, and translational invariance,
      !! which is sensitive to a different set of mistakes and costs nothing.
      !!
      !! The tolerance is loose and cannot be tightened. The gradient being
      !! differenced carries the grid-response terms that `xc_hessian` omits, as
      !! PySCF's does by default, so the two disagree by that term however exact
      !! the rest is. It is far above step noise and far below what an assembly
      !! error produces -- a missing kernel in the relaxed mean field was 0.26
      !! on this system, and a missing potential derivative 0.87.
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp) :: shifted(3, 3), fd, worst, rowsum, wtrans
      real(dp), allocatable :: hess(:, :, :, :), plus(:, :), minus(:, :)
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b

      if (.not. xc_available()) return

      call ks_end_to_end_for("lda_x", error)
      if (allocated(error)) return
      call ks_end_to_end_for("pbe", error)
      if (allocated(error)) return
      call ks_end_to_end_for("b3lyp", error)
      if (allocated(error)) return
      call ks_end_to_end_for("tpss", error)
      if (allocated(error)) return
      ! Named in the coverage table and, until this line, not exercised by
      ! anything.
      call ks_end_to_end_for("m06-l", error)
      if (allocated(error)) return
      ! The `-V` case exercises all three VV10 pieces through the assembled
      ! Hessian: the explicit second derivative, the perturbed Fock and the
      ! response kernel. The gradient being differenced carries the NLC
      ! grid's response as well as the semilocal grid's, so the residual is
      ! larger than the pieces' own step errors for the reason the header
      ! sets out -- measured at 1.2e-3 here (invariance 1.0e-3), the size of
      ! the TPSS-class cases and inside the shared tolerance.
      call ks_end_to_end_for("b97m-v", error)
   end subroutine test_ks_end_to_end

   subroutine ks_end_to_end_for(functional, error)
      !! One functional, assembled Hessian against a difference of its gradient
      character(len=*), intent(in) :: functional
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: H = 1.0e-4_dp
      real(dp) :: shifted(3, 3), fd, worst, rowsum, wtrans
      real(dp), allocatable :: hess(:, :, :, :), plus(:, :), minus(:, :)
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b

      call ks_at_f(WATER, functional, hess, err, ok)
      call check(error, ok, "the Kohn-Sham Hessian failed for "//functional)
      if (allocated(error)) return

      wtrans = 0.0_dp
      do a = 1, 3
         do b = 1, 3
            do ia = 1, 3
               rowsum = sum(hess(a, b, ia, :))
               wtrans = max(wtrans, abs(rowsum))
            end do
         end do
      end do

      worst = 0.0_dp
      do ia = 1, 3
         do a = 1, 3
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) + H
            call dft_gradient_at(shifted, functional, plus, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            shifted = WATER
            shifted(a, ia) = shifted(a, ia) - H
            call dft_gradient_at(shifted, functional, minus, err, ok)
            if (.not. ok) then
               call check(error, .false., "a displaced gradient failed")
               return
            end if
            do ja = 1, 3
               do b = 1, 3
                  fd = (plus(b, ja) - minus(b, ja))/(2.0_dp*H)
                  worst = max(worst, abs(hess(a, b, ia, ja) - fd))
               end do
            end do
         end do
      end do
      ! Loose because the omitted grid response breaks it by an amount that
      ! depends on the functional and the grid, not because anything here is
      ! approximate. On this system PySCF's own Hessian violates it by 5.4e-4
      ! at LDA, 5.7e-4 at PBE, 1.2e-3 at TPSS and 6.7e-3 at M06-L, and ours
      ! agrees with PySCF's to 1e-8 throughout. A gross assembly error shows up
      ! here at 1e-1 and larger, which is what this is guarding against.
      call check(error, wtrans < 2.0e-2_dp, &
                 "the Kohn-Sham Hessian is not translationally invariant")
      if (allocated(error)) return
      ! Bounded by the finite difference, not by the Hessian. PySCF's own
      ! analytic Kohn-Sham Hessian is 4.3e-4 from a difference of PySCF's own
      ! gradient on this system, and ours is 5.2e-4 from a difference of ours --
      ! the same quadrature noise, not an error in either. The tight check is
      ! against PySCF's analytic Hessian directly, where this agrees to 1.8e-8;
      ! that comparison lives in the validation suite because it needs PySCF.
      call check(error, worst < 1.0e-2_dp, &
                 "the Kohn-Sham Hessian disagrees with differences of the gradient for "// &
                 functional)
   end subroutine ks_end_to_end_for

   subroutine ks_at_f(coords, functional, hess, err, ok)
      real(dp), intent(in) :: coords(3, 3)
      character(len=*), intent(in) :: functional
      real(dp), allocatable, intent(out) :: hess(:, :, :, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok
      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, functional, ctx, err, level=3, allow_half=.true.)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (.not. err%has_error()) then
         call ks_hessian(mol, WATER_Z, scf%density, scf%orbitals, scf%orbital_energies, &
                         5, ctx, ctx%exx_fraction, hess, err)
      end if
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine ks_at_f

   subroutine dft_gradient_at(coords, functional, gradient, err, ok)
      real(dp), intent(in) :: coords(3, 3)
      character(len=*), intent(in) :: functional
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, functional, ctx, err, level=3, allow_half=.true.)
      if (err%has_error()) return
      call run_libcint_rhf(mol, 10, 100, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err, xc=ctx)
      if (.not. err%has_error()) then
         call libcint_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                                   orbital_energies=scf%orbital_energies, &
                                   n_occupied=5, gradient=gradient, error=err, xc=ctx)
      end if
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine dft_gradient_at

   subroutine test_kernel_deriv(error)
      !! The kernel's nuclear derivative, against central differences of
      !! `xc_kernel_apply` across geometry
      !!
      !! `xc_kernel_deriv` is `d/dR` of the kernel applied to a trial density,
      !! with both matrices and the grid held fixed -- so differencing
      !! `xc_kernel_apply` itself over displaced nuclei, on the one context
      !! and the one pair of matrices, is exact up to step error. Every entry
      !! of every `(3, natm)` component is compared; an aggregate would pass
      !! on errors that cancel across components.
      !!
      !! Both families, because an LDA-only implementation gets the rho
      !! channel right while every sigma channel stays silently zero: B88
      !! exercises all of `g_xc`, `lda_x` only `v3rho3`.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst

      if (.not. xc_available()) return

      call kernel_deriv_case(WATER_Z, WATER_SYM, WATER, "sto-3g", "gga_x_b88", &
                             1.0e-4_dp, worst, error)
      if (allocated(error)) return
      call check(error, worst < 1.0e-10_dp, &
                 "gga_x_b88: the kernel derivative disagrees with differences "// &
                 "of the kernel apply")
         !! Eight times the measured worst, 1.27e-11 at H = 1e-4. That is step
         !! error and nothing else: 3.16e-12 at 5e-5 and 5.08e-11 at 2e-4,
         !! ratios 4.02 and 4.00. A break test zeroing the `g_xc` group misses
         !! by 2.08e-2 at every step -- nine orders up and step-independent.
      if (allocated(error)) return

      call kernel_deriv_case(WATER_Z, WATER_SYM, WATER, "sto-3g", "lda_x", &
                             1.0e-4_dp, worst, error)
      if (allocated(error)) return
      call check(error, worst < 1.0e-10_dp, &
                 "lda_x: the kernel derivative disagrees with differences "// &
                 "of the kernel apply")
         !! Eight times the measured worst, 1.30e-11 at H = 1e-4, scaling as
         !! H^2 (3.24e-12 at 5e-5 and 5.21e-11 at 2e-4, ratios 4.02 and 4.00).
         !! The same break test misses by 1.95e-2 here, so the rho channel
         !! alone pins `v3rho3` too.
   end subroutine test_kernel_deriv

   subroutine test_kernel_deriv_dimer(error)
      !! The same rung on the water dimer in 6-31G(d): d functions, six atoms,
      !! and no reflection plane for a wrong component to hide behind
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst

      if (.not. xc_available()) return

      call kernel_deriv_case(DIMER_Z6, DIMER_SYM6, DIMER_XYZ6, "6-31g*", &
                             "gga_x_b88", 1.0e-4_dp, worst, error)
      if (allocated(error)) return
      call check(error, worst < 2.0e-9_dp, &
                 "the kernel derivative disagrees with differences of the "// &
                 "kernel apply on a dimer")
         !! Eight times the measured worst, 2.62e-10 at H = 1e-4 across all
         !! 18 x nao x nao entries; 6.54e-11 at 5e-5, ratio 4.01. The `g_xc`
         !! break test misses by 8.31e-2 here, step-independent.
   end subroutine test_kernel_deriv_dimer

   subroutine kernel_deriv_case(z, sym, coords, basis, functional, hstep, worst, error)
      !! One geometry and functional: every entry of `xc_kernel_deriv` against
      !! a central difference of `xc_kernel_apply` over displaced nuclei
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, functional
      real(dp), intent(in) :: hstep
      real(dp), intent(out) :: worst
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), ta(:, :), tb(:, :)
      real(dp), allocatable :: h1(:, :, :, :), vp(:, :), vm(:, :), shifted(:, :)
      logical :: ok
      integer :: natm, nao, ia, a, u, v

      worst = huge(1.0_dp)
      natm = size(z)

      call scf_for(z, sym, coords, basis, functional, ctx, dens, err, ok)
      call check(error, ok, "the "//functional//" reference failed")
      if (allocated(error)) return
      nao = size(dens, 1)
      call kernel_trials(dens, ta, tb)

      allocate (h1(nao, nao, 3, natm))
      h1 = 0.0_dp
      call build_libcint_molecule(z, sym, coords, basis, mol, err)
      if (.not. err%has_error()) call xc_kernel_deriv(ctx, mol, dens, ta, h1, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the kernel derivative failed")
      if (allocated(error)) return

      allocate (shifted(3, natm))
      worst = 0.0_dp
      do ia = 1, natm
         do a = 1, 3
            shifted = coords
            shifted(a, ia) = coords(a, ia) + hstep
            call kernel_matrix_at(ctx, z, sym, shifted, basis, dens, ta, vp, err)
            shifted(a, ia) = coords(a, ia) - hstep
            call kernel_matrix_at(ctx, z, sym, shifted, basis, dens, ta, vm, err)
            if (err%has_error()) exit
            do v = 1, nao
               do u = 1, nao
                  worst = max(worst, abs(h1(u, v, a, ia) &
                                         - (vp(u, v) - vm(u, v))/(2.0_dp*hstep)))
               end do
            end do
         end do
         if (err%has_error()) exit
      end do
      call check(error,.not. err%has_error(), "a displaced kernel apply failed")
   end subroutine kernel_deriv_case

   subroutine test_kernel2(error)
      !! The second kernel, against central differences of `xc_kernel_apply`
      !! in the *density*
      !!
      !! `xc_kernel2_apply(D, ta, tb)` is the derivative of
      !! `xc_kernel_apply(D, tb)` when `D` moves by `ta`, at fixed geometry --
      !! so `(V[D + h ta](tb) - V[D - h ta](tb)) / 2h` with nothing else
      !! moving is exact up to step error, entry by entry. The two trials are
      !! deliberately unlike each other; a and b enter symmetrically, so a
      !! transposed channel pairing has to disagree with an asymmetric pair.
      !!
      !! Both families, for the same reason as the derivative test above: an
      !! LDA-only second kernel is right by accident on `lda_x` and wrong on
      !! every sigma channel of B88.
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst

      if (.not. xc_available()) return

      ! The step is 1e-3 where the geometric tests use 1e-4, because this
      ! comparison is far more converged: at 1e-4 the disagreement bottoms out
      ! at the roundoff floor, ~1e-13, and stops scaling -- there is no step
      ! error left to measure down there.
      call kernel2_case(WATER_Z, WATER_SYM, WATER, "sto-3g", "gga_x_b88", &
                        1.0e-3_dp, worst, error)
      if (allocated(error)) return
      call check(error, worst < 1.0e-10_dp, &
                 "gga_x_b88: the second kernel disagrees with differences "// &
                 "of the kernel apply in the density")
         !! Eight times the measured worst, 1.16e-11 at H = 1e-3, which is
         !! step error and nothing else: 2.90e-12 at 5e-4 and 4.66e-11 at
         !! 2e-3, ratios 4.02 and 4.00. A break test zeroing the `g_xc` group
         !! misses by 5.67e-3 at every step -- eight orders up and
         !! step-independent.
      if (allocated(error)) return

      call kernel2_case(WATER_Z, WATER_SYM, WATER, "sto-3g", "lda_x", &
                        1.0e-3_dp, worst, error)
      if (allocated(error)) return
      call check(error, worst < 1.0e-10_dp, &
                 "lda_x: the second kernel disagrees with differences "// &
                 "of the kernel apply in the density")
         !! Eight times the measured worst, 1.18e-11 at H = 1e-3, scaling as
         !! H^2 (2.96e-12 at 5e-4 and 4.73e-11 at 2e-3, ratios 4.00 and 4.00).
         !! The same break test misses by 4.22e-3 here: for an LDA the second
         !! kernel *is* the `g_xc` term, so dropping it returns zero.
   end subroutine test_kernel2

   subroutine test_kernel2_dimer(error)
      !! The second kernel on the water dimer in 6-31G(d), where d functions
      !! and the absence of any symmetry plane leave an error nowhere to hide
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: worst

      if (.not. xc_available()) return

      call kernel2_case(DIMER_Z6, DIMER_SYM6, DIMER_XYZ6, "6-31g*", &
                        "gga_x_b88", 1.0e-4_dp, worst, error)
      if (allocated(error)) return
      call check(error, worst < 1.5e-9_dp, &
                 "the second kernel disagrees with differences of the kernel "// &
                 "apply in the density on a dimer")
         !! Seven times the measured worst, 2.21e-10 at H = 1e-4 across every
         !! entry; 5.4-5.7e-11 at 5e-5 over repeated runs, ratio 3.9-4.1 --
         !! the smaller step carries a few per cent of run-to-run scatter from
         !! the threaded kernel's merge order, the 1e-4 figure repeats to four
         !! digits. The `g_xc` break test misses by 1.16e-1 here,
         !! step-independent.
   end subroutine test_kernel2_dimer

   subroutine kernel2_case(z, sym, coords, basis, functional, hstep, worst, error)
      !! One geometry and functional: every entry of `xc_kernel2_apply`
      !! against `(V[D + h ta](tb) - V[D - h ta](tb)) / 2h` at fixed geometry
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, functional
      real(dp), intent(in) :: hstep
      real(dp), intent(out) :: worst
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(xc_context_t) :: ctx
      type(error_t) :: err
      real(dp), allocatable :: dens(:, :), ta(:, :), tb(:, :)
      real(dp), allocatable :: v2(:, :), vp(:, :), vm(:, :)
      logical :: ok
      integer :: nao, u, v

      worst = huge(1.0_dp)

      call scf_for(z, sym, coords, basis, functional, ctx, dens, err, ok)
      call check(error, ok, "the "//functional//" reference failed")
      if (allocated(error)) return
      nao = size(dens, 1)
      call kernel_trials(dens, ta, tb)

      allocate (v2(nao, nao))
      v2 = 0.0_dp
      call build_libcint_molecule(z, sym, coords, basis, mol, err)
      if (.not. err%has_error()) call xc_kernel2_apply(ctx, mol, dens, ta, tb, v2, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the second kernel apply failed")
      if (allocated(error)) return

      call kernel_matrix_at(ctx, z, sym, coords, basis, dens + hstep*ta, tb, vp, err)
      call kernel_matrix_at(ctx, z, sym, coords, basis, dens - hstep*ta, tb, vm, err)
      call check(error,.not. err%has_error(), "a displaced-density kernel apply failed")
      if (allocated(error)) return

      worst = 0.0_dp
      do v = 1, nao
         do u = 1, nao
            worst = max(worst, abs(v2(u, v) - (vp(u, v) - vm(u, v))/(2.0_dp*hstep)))
         end do
      end do
   end subroutine kernel2_case

   subroutine kernel_matrix_at(ctx, z, sym, coords, basis, dens, dtilde, v, err)
      !! `xc_kernel_apply` into a zeroed matrix: the quantity both kernel
      !! tests difference, at whatever geometry and density they hand in --
      !! the grid never moves because the context is passed in, not rebuilt
      type(xc_context_t), intent(inout) :: ctx
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis
      real(dp), intent(in) :: dens(:, :), dtilde(:, :)
      real(dp), allocatable, intent(out) :: v(:, :)
      type(error_t), intent(inout) :: err

      type(libcint_molecule_t) :: mol

      if (err%has_error()) return
      allocate (v(size(dens, 1), size(dens, 2)))
      v = 0.0_dp
      call build_libcint_molecule(z, sym, coords, basis, mol, err)
      if (err%has_error()) return
      call xc_kernel_apply(ctx, mol, dens, dtilde, v, err)
      call mol%destroy()
   end subroutine kernel_matrix_at

   subroutine kernel_trials(dens, ta, tb)
      !! Two symmetric trial densities, deliberately unlike each other and
      !! neither proportional to the density -- a proportional trial moves rho
      !! and sigma together and cannot separate the channels, and identical
      !! trials would let a transposed a/b pairing pass unseen
      real(dp), intent(in) :: dens(:, :)
      real(dp), allocatable, intent(out) :: ta(:, :), tb(:, :)

      integer :: n, i, j

      n = size(dens, 1)
      allocate (ta(n, n), tb(n, n))
      do j = 1, n
         do i = 1, n
            ta(i, j) = 0.05_dp/(1.0_dp + real(abs(i - j), dp)) + 0.01_dp*dens(i, j)
            tb(i, j) = 0.1_dp*cos(real(i - j, dp)) + 0.03_dp*sin(real(i + j, dp))
         end do
         tb(j, j) = tb(j, j) + 0.2_dp
      end do
      ta = 0.5_dp*(ta + transpose(ta))
      tb = 0.5_dp*(tb + transpose(tb))
   end subroutine kernel_trials

   subroutine scf_for(z, sym, coords, basis, functional, ctx, dens, err, ok)
      !! One converged Kohn-Sham reference on an arbitrary geometry: the grid
      !! and the density the kernel tests hold fixed
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: coords(:, :)
      character(len=*), intent(in) :: basis, functional
      type(xc_context_t), intent(out) :: ctx
      real(dp), allocatable, intent(out) :: dens(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      ok = .false.
      call build_libcint_molecule(z, sym, coords, basis, mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, functional, ctx, err, level=3, allow_half=.true.)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if
      call run_libcint_rhf(mol, sum(z), 100, 1.0e-12_dp, 1.0e-10_dp, .false., &
                           scf, err, xc=ctx)
      if (.not. err%has_error()) dens = scf%density
      call mol%destroy()
      ok = .not. err%has_error()
   end subroutine scf_for

end module test_mqc_xc_hessian

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_xc_hessian, only: collect_mqc_xc_hessian
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_xc_hessian", collect_mqc_xc_hessian)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
