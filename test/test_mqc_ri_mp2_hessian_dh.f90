!! A double hybrid's fitted perturbative Hessian
module test_mqc_ri_mp2_hessian_dh
   !! `ri_mp2_correlation_hessian` with a functional is the correlation block of a
   !! double hybrid over a Kohn-Sham reference, scaled by its PT2 coefficient. Two
   !! checks, each against something that owes it nothing:
   !!
   !! * **Its own gradient.** One column of the block against a seven-point
   !!   stencil of `czt_ri_mp2_gradient` with the same functional, on one grid
   !!   held fixed -- the grid the analytic Hessian holds fixed -- and with an
   !!   exact reference operator, `fit_response = .false.`, so the energy being
   !!   differentiated is the one the Hessian differentiates.
   !! * **The conventional block.** The same double hybrid through
   !!   `mp2_correlation_hessian`. The two differ by the fitting error and by
   !!   nothing else, so the relative difference has to be small and, as a guard
   !!   against the fitted path quietly being the exact one, not zero. The
   !!   dipole derivatives are held to the same standard.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_xc, only: xc_context_t, xc_context_create, xc_available
   use mqc_czt_ri_mp2_gradient, only: czt_ri_mp2_gradient
   use mqc_czt_ri_mp2_hessian, only: ri_mp2_correlation_hessian
   use mqc_czt_mp2_hessian, only: mp2_correlation_hessian
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_ri_mp2_hessian_dh_tests

   integer, parameter :: Z(3) = [8, 1, 1]
   character(len=2), parameter :: SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.05_dp, 0.02_dp, -0.03_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.04_dp, 1.4308_dp, 1.1078_dp], [3, 3])
   character(len=*), parameter :: FUNCTIONAL = "b2plyp"
   real(dp), parameter :: STEP = 2.0e-3_dp
   real(dp), parameter :: TOL_FD = 5.0e-8_dp
      !! Measured 2e-9 on one thread; a missing or mis-signed piece of the kernel
      !! response is orders larger.
   real(dp), parameter :: FITTING_LIMIT = 5.0e-3_dp
      !! Relative, RI against conventional. Measured 1e-3 or less with
      !! cc-pVDZ-RIFIT on water; a wrong block is order one.
   integer, parameter :: PERT_ATOM = 1
   integer, parameter :: PERT_CART = 3

contains

   subroutine collect_mqc_ri_mp2_hessian_dh_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("b2plyp_against_the_fitted_gradient", against_gradient), &
                  new_unittest("b2plyp_6_31g_against_the_conventional_block", &
                               against_conventional_631g), &
                  new_unittest("b2plyp_cc_pvdz_against_the_conventional_block", &
                               against_conventional_dz), &
                  new_unittest("batching_changes_nothing", batching) &
                  ]
   end subroutine collect_mqc_ri_mp2_hessian_dh_tests

   subroutine setup_case(basis, aux_basis, geo, mol, aux, ctx, scf, err)
      character(len=*), intent(in) :: basis, aux_basis
      real(dp), intent(in) :: geo(3, 3)
      type(czt_molecule_t), intent(out) :: mol, aux
      type(xc_context_t), intent(inout) :: ctx
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      call build_czt_molecule(Z, SYM, geo, basis, mol, err)
      if (err%has_error()) return
      call build_czt_molecule(Z, SYM, geo, aux_basis, aux, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 10, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err, xc=ctx)
   end subroutine setup_case

   subroutine make_context(basis, ctx, err)
      character(len=*), intent(in) :: basis
      type(xc_context_t), intent(out) :: ctx
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol

      call build_czt_molecule(Z, SYM, WATER, basis, mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, FUNCTIONAL, ctx, err, level=3)
      call mol%destroy()
   end subroutine make_context

   subroutine fitted_pt2_gradient(ctx, basis, aux_basis, coords, gradient, err)
      type(xc_context_t), intent(inout) :: ctx
      character(len=*), intent(in) :: basis, aux_basis
      real(dp), intent(in) :: coords(3, 3)
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf

      call setup_case(basis, aux_basis, coords, mol, aux, ctx, scf, err)
      if (err%has_error()) return
      call czt_ri_mp2_gradient(mol, aux, scf%orbitals, scf%orbital_energies, 5, gradient, err, &
                               xc=ctx, scf_density=scf%density, pt2_scale=ctx%pt2_fraction, &
                               fit_response=.false., fixed_grid=.true.)
      call mol%destroy()
      call aux%destroy()
   end subroutine fitted_pt2_gradient

   subroutine against_gradient(error)
      !! One column of the block, `d(gradient)/d(atom 1, z)`, by the seven-point
      !! stencil
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: WEIGHT(6) = [-1.0_dp, 9.0_dp, -45.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]
      integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
      type(czt_molecule_t) :: mol, aux
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: hess(:, :, :, :), gradient(:, :)
      real(dp) :: coords(3, 3), acc(3, 3), fd(9), an(9)
      integer :: point, ia, a, threads

      if (.not. xc_available()) return
      ! One thread, as the conventional gate does: the stencil divides by 60h and the
      ! unordered OpenMP merges do not sum in a fixed order.
      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call make_context("6-31g", ctx, err)
      call check(error, .not. err%has_error(), "the functional context must build")
      if (allocated(error)) return

      acc = 0.0_dp
      do point = 1, 6
         coords = WATER
         coords(PERT_CART, PERT_ATOM) = coords(PERT_CART, PERT_ATOM) + real(OFFSET(point), dp)*STEP
         call fitted_pt2_gradient(ctx, "6-31g", "cc-pvdz-rifit", coords, gradient, err)
         call check(error, .not. err%has_error(), "a displaced gradient failed: "// &
                    err%get_full_trace())
         if (allocated(error)) then
            call omp_set_num_threads(threads)
            return
         end if
         acc = acc + WEIGHT(point)*gradient
      end do
      fd = reshape(acc/(60.0_dp*STEP), [9])

      call setup_case("6-31g", "cc-pvdz-rifit", WATER, mol, aux, ctx, scf, err)
      if (.not. err%has_error()) &
         call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, 5, hess, &
                                         err, tol=1.0e-13_dp, xc=ctx, scf_density=scf%density, &
                                         pt2_scale=ctx%pt2_fraction)
      call omp_set_num_threads(threads)
      call check(error, .not. err%has_error(), "the Hessian must run: "//err%get_full_trace())
      if (allocated(error)) return
      do ia = 1, 3
         do a = 1, 3
            an(3*(ia - 1) + a) = hess(a, PERT_CART, ia, PERT_ATOM)
         end do
      end do
      write (*, "(a,es10.3)") "        max |analytic - fd| = ", maxval(abs(an - fd))
      call check(error, maxval(abs(an - fd)) < TOL_FD, &
                 "the fitted double-hybrid Hessian does not difference its own gradient")
      call mol%destroy()
      call aux%destroy()
   end subroutine against_gradient

   subroutine conventional_case(basis, aux_basis, error)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: basis, aux_basis

      type(czt_molecule_t) :: mol, aux
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: ri(:, :, :, :), conv(:, :, :, :), discard(:, :, :, :)
      real(dp), allocatable :: dip_ri(:, :), dip_conv(:, :)
      real(dp) :: rel, rel_dip

      if (.not. xc_available()) return
      call make_context(basis, ctx, err)
      if (.not. err%has_error()) call setup_case(basis, aux_basis, WATER, mol, aux, ctx, scf, err)
      call check(error, .not. err%has_error(), "the Kohn-Sham reference must converge: "// &
                 err%get_full_trace())
      if (allocated(error)) return

      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, 5, ri, err, &
                                      tol=1.0e-13_dp, xc=ctx, scf_density=scf%density, &
                                      pt2_scale=ctx%pt2_fraction, dipole_derivatives=dip_ri)
      if (.not. err%has_error()) &
         call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, scf%density, 5, 0, &
                                      conv, discard, err, xc=ctx, pt2_scale=ctx%pt2_fraction, &
                                      dipole_derivatives=dip_conv)
      call check(error, .not. err%has_error(), "both blocks must run: "//err%get_full_trace())
      if (allocated(error)) return

      rel = maxval(abs(ri - conv))/maxval(abs(conv))
      rel_dip = maxval(abs(dip_ri - dip_conv))/maxval(abs(dip_conv))
      write (*, "(a,a,2es10.3)") "        ", basis//": relative RI vs conventional, block and dipole ", &
         rel, rel_dip
      call check(error, rel < FITTING_LIMIT, "the fitted block is not the conventional one "// &
                 "to the fitting error")
      if (allocated(error)) return
      call check(error, rel > 1.0e-7_dp, "the fitted block equals the exact one: it did not fit")
      if (allocated(error)) return
      call check(error, rel_dip < FITTING_LIMIT, "the fitted dipole derivative is not the "// &
                 "conventional one to the fitting error")
      call mol%destroy()
      call aux%destroy()
   end subroutine conventional_case

   subroutine against_conventional_631g(error)
      type(error_type), allocatable, intent(out) :: error

      call conventional_case("6-31g", "cc-pvdz-rifit", error)
   end subroutine against_conventional_631g

   subroutine against_conventional_dz(error)
      type(error_type), allocatable, intent(out) :: error

      call conventional_case("cc-pvdz", "cc-pvdz-rifit", error)
   end subroutine against_conventional_dz

   subroutine batching(error)
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(xc_context_t) :: ctx
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: whole(:, :, :, :), split(:, :, :, :)

      if (.not. xc_available()) return
      call make_context("6-31g", ctx, err)
      if (.not. err%has_error()) &
         call setup_case("6-31g", "cc-pvdz-rifit", WATER, mol, aux, ctx, scf, err)
      call check(error, .not. err%has_error(), "the reference must converge")
      if (allocated(error)) return
      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, 5, whole, &
                                      err, tol=1.0e-13_dp, xc=ctx, scf_density=scf%density, &
                                      pt2_scale=ctx%pt2_fraction)
      call ri_mp2_correlation_hessian(mol, aux, scf%orbitals, scf%orbital_energies, 5, split, &
                                      err, tol=1.0e-13_dp, xc=ctx, scf_density=scf%density, &
                                      pt2_scale=ctx%pt2_fraction, batch=2, atom_batch=1)
      call check(error, .not. err%has_error(), "both runs must succeed: "//err%get_full_trace())
      if (allocated(error)) return
      call check(error, maxval(abs(whole - split)) < 1.0e-10_dp, "the batch sizes change the block")
      call mol%destroy()
      call aux%destroy()
   end subroutine batching

end module test_mqc_ri_mp2_hessian_dh

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ri_mp2_hessian_dh, only: collect_mqc_ri_mp2_hessian_dh_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_ri_mp2_hessian_dh", collect_mqc_ri_mp2_hessian_dh_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
