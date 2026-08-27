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
                             vv10_kernel_apply
   use mqc_libcint_hessian, only: ks_hessian
   use mqc_libcint_gradient, only: libcint_scf_gradient, vv10_gradient_fixed_grid
   use mqc_libcint_xc_hessian, only: xc_hessian, xc_gradient_fixed_grid, &
                                     xc_potential_deriv, vv10_hessian, &
                                     vv10_potential_deriv
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

contains

   subroutine collect_mqc_xc_hessian(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("xc_hessian_differences_its_own_gradient", test_against_fd), &
                  new_unittest("xc_hessian_mgga_differences_its_gradient", test_gga_against_fd), &
                  new_unittest("vv10_fixed_grid_gradient_differences_the_energy", &
                               test_vv10_fixed_grid), &
                  new_unittest("vv10_hessian_differences_the_fixed_grid_gradient", &
                               test_vv10_hessian), &
                  new_unittest("xc_potential_derivative_matches_differences", test_vxc_deriv), &
                  new_unittest("vv10_potential_derivative_matches_differences", &
                               test_vv10_vxc_deriv), &
                  new_unittest("vv10_kernel_differences_the_potential_in_density", &
                               test_vv10_kernel_apply), &
                  new_unittest("ks_hessian_differences_the_dft_gradient", test_ks_end_to_end) &
                  ]
   end subroutine collect_mqc_xc_hessian

   subroutine xc_at(ctx, coords, gradient, hess, density, err, ok)
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
      real(dp), intent(out), optional :: gradient(3, 3)
      real(dp), intent(out), optional :: hess(3, 3, 3, 3)
      real(dp), intent(in) :: density(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, "sto-3g", mol, err)
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

   subroutine reference_state(functional, ctx, density, err, ok)
      !! One converged reference: the grid and the density everything below uses
      character(len=*), intent(in) :: functional
      type(xc_context_t), intent(out) :: ctx
      real(dp), allocatable, intent(out) :: density(:, :)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, functional, ctx, err, level=3)
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

      real(dp), parameter :: H = 1.0e-3_dp
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

      real(dp), parameter :: H = 1.0e-3_dp
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

      real(dp), parameter :: H = 2.5e-4_dp
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

      real(dp), parameter :: H = 1.0e-3_dp
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

      real(dp), parameter :: H = 1.0e-3_dp
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

      real(dp), parameter :: H = 1.0e-3_dp
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

      real(dp), parameter :: H = 1.0e-3_dp
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

      real(dp), parameter :: H = 1.0e-3_dp
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

      real(dp), parameter :: H = 1.0e-3_dp
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
      call xc_context_create(mol, functional, ctx, err, level=3)
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
      call xc_context_create(mol, functional, ctx, err, level=3)
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
