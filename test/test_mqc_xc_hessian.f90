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
                             xc_add_potential
   use mqc_libcint_xc_hessian, only: xc_hessian, xc_gradient_fixed_grid, &
                                     xc_potential_deriv
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
                  new_unittest("xc_hessian_gga_differences_its_gradient", test_gga_against_fd), &
                  new_unittest("xc_hessian_refuses_a_meta_gga", test_refuses_mgga), &
                  new_unittest("xc_potential_derivative_matches_differences", test_vxc_deriv) &
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

      real(dp), parameter :: H = 2.0e-3_dp
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

      real(dp), parameter :: H = 2.0e-3_dp
      real(dp), parameter :: TOL = 2.0e-5_dp
      real(dp) :: hess(3, 3, 3, 3), plus(3, 3), minus(3, 3), shifted(3, 3)
      real(dp) :: fd
      real(dp), allocatable :: dens(:, :)
      type(xc_context_t) :: ctx
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, ja, b

      if (.not. xc_available()) return

      call reference_state("gga_x_pbe", ctx, dens, err, ok)
      call check(error, ok, "the reference Kohn-Sham state failed")
      if (allocated(error)) return

      call xc_at(ctx, WATER, hess=hess, density=dens, err=err, ok=ok)
      call check(error, ok, "the analytic GGA exchange-correlation Hessian failed")
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
                             more="GGA exchange-correlation Hessian entry disagrees "// &
                             "with a difference of its own gradient")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do
   end subroutine test_gga_against_fd

   subroutine test_refuses_mgga(error)
      !! A meta-GGA is refused rather than silently missing its tau terms
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: hess(3, 3, 3, 3)
      real(dp), allocatable :: dens(:, :)
      type(xc_context_t) :: ctx
      type(error_t) :: err
      logical :: ok

      if (.not. xc_available()) return

      call reference_state("mgga_x_scan", ctx, dens, err, ok)
      if (.not. ok) return
      call xc_at(ctx, WATER, hess=hess, density=dens, err=err, ok=ok)
      call check(error,.not. ok, "a meta-GGA Hessian must be refused, not approximated")
   end subroutine test_refuses_mgga

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

      real(dp), parameter :: H = 2.0e-3_dp
      real(dp) :: shifted(3, 3), worst
      real(dp), allocatable :: dens(:, :), h1(:, :, :, :), vp(:, :), vm(:, :)
      type(xc_context_t) :: ctx
      type(libcint_molecule_t) :: mol
      type(error_t) :: err
      logical :: ok
      integer :: ia, a, nao, u, v

      if (.not. xc_available()) return

      call reference_state("lda_x", ctx, dens, err, ok)
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
                 "differences of the potential")
   end subroutine test_vxc_deriv

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
