module test_mqc_czt_df_hessian
   !! The analytic Hessian of a density-fitted Hartree-Fock energy.
   !!
   !! Held to the one thing it has to be: the derivative of the fitted
   !! gradient. A five-point stencil of `czt_scf_gradient` with the same
   !! auxiliary basis is a reference that shares no line with the Hessian --
   !! the gradient contracts first-derivative integrals against densities of
   !! its own, and knows nothing of the coupled-perturbed solve, the
   !! second-derivative integrals or the product term this adds.
   !!
   !! That matters most for the term a fitted Hessian can lose without looking
   !! wrong: the response of the fit itself to a moved auxiliary basis. Carried
   !! only to first order -- PySCF's `auxbasis_response = 1` -- the Hessian is
   !! still symmetric and still translationally invariant, and is 1e-4 of an
   !! element from the differenced gradient on water. The validation suite holds the
   !! same Hessian to PySCF's; this holds it to our own gradient, which is what
   !! a geometry optimisation and a differenced Hessian actually use.
   !!
   !! The geometry is off every axis and plane, so no block is zero for a
   !! reason that has nothing to do with the integrals.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_gradient, only: czt_scf_gradient
   use mqc_czt_hessian, only: df_rhf_hessian, rhf_hessian, response_controls_t
   use mqc_czt_df_hessian, only: df_metric_dropped
   implicit none
   private

   public :: collect_mqc_czt_df_hessian_tests

   integer, parameter :: N_ATOMS = 3
   integer, parameter :: N_COORD = 3*N_ATOMS
   integer, parameter :: N_OCC = 5
   integer, parameter :: N_ELECTRONS = 10
   real(dp), parameter :: WATER(3, N_ATOMS) = reshape( &
                          [0.05_dp, 0.02_dp, -0.03_dp, &
                           0.0_dp, -1.4308_dp, 1.1078_dp, &
                           0.04_dp, 1.4308_dp, 1.1078_dp], [3, N_ATOMS])
   character(len=2), parameter :: SYM(N_ATOMS) = ["O ", "H ", "H "]
   integer, parameter :: NUMBERS(N_ATOMS) = [8, 1, 1]

   !! Bohr. Small enough for the stencil's own error to sit near 1e-9, large
   !! enough that the SCF's convergence does not show through the division.
   real(dp), parameter :: STEP = 2.0e-3_dp
   !! Hartree/Bohr^2. The stencil, not the Hessian, sets this.
   real(dp), parameter :: TOL_FD = 3.0e-7_dp
   !! What separates a fitted Hessian from an exact one on this molecule, at
   !! the least: the two must differ, or the fit was never differentiated.
   real(dp), parameter :: MIN_FIT_EFFECT = 1.0e-5_dp
   !! On the SCF's commutator, so the density and not only the energy is
   !! converged.
   real(dp), parameter :: SCF_GRAD_TOL = 1.0e-10_dp

contains

   subroutine collect_mqc_czt_df_hessian_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("df_hessian_is_the_fitted_gradients_derivative", &
                               against_gradient_differences), &
                  new_unittest("df_hessian_with_d_functions", &
                               against_gradient_differences_polarised), &
                  new_unittest("df_hessian_is_not_the_exact_one", differs_from_exact), &
                  new_unittest("df_hessian_metric_keeps_every_mode", metric_is_whole) &
                  ]
   end subroutine collect_mqc_czt_df_hessian_tests

   subroutine converge(geo, basis, aux_basis, mol, aux, scf, err)
      !! A fitted closed-shell SCF on water at `geo`, converged hard
      real(dp), intent(in) :: geo(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      type(czt_molecule_t), intent(out) :: mol, aux
      type(rhf_result_t), intent(out) :: scf
      type(error_t), intent(inout) :: err

      call build_czt_molecule(NUMBERS, SYM, geo, basis, mol, err)
      if (err%has_error()) return
      call build_czt_molecule(NUMBERS, SYM, geo, aux_basis, aux, err)
      if (err%has_error()) return
      ! The commutator is named because the density is what gets used: left to
      ! derive itself from the energy tolerance it stops near 1e-7, and the
      ! stencil below then divides the difference between two such stops by a
      ! step of 2e-3. Without it this test passed on one thread and failed on
      ! four about one run in two, by up to 1e-6.
      call run_czt_rhf(mol, N_ELECTRONS, 300, 1.0e-14_dp, 1.0e-12_dp, .false., scf, err, &
                       aux=aux, grad_tol=SCF_GRAD_TOL)
   end subroutine converge

   subroutine fitted_gradient(geo, basis, aux_basis, grad, err)
      !! The fitted SCF gradient at one geometry, flattened
      real(dp), intent(in) :: geo(:, :)
      character(len=*), intent(in) :: basis, aux_basis
      real(dp), intent(out) :: grad(:)
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      real(dp), allocatable :: g(:, :)

      call converge(geo, basis, aux_basis, mol, aux, scf, err)
      if (err%has_error()) return
      call czt_scf_gradient(mol, scf%density, orbitals=scf%orbitals, &
                            orbital_energies=scf%orbital_energies, n_occupied=N_OCC, &
                            gradient=g, error=err, aux=aux)
      if (err%has_error()) return
      grad = reshape(g, [size(grad)])
      call mol%destroy()
      call aux%destroy()
   end subroutine fitted_gradient

   subroutine to_matrix(h4, matrix)
      !! `(3, 3, natm, natm)` to `(3N, 3N)`, row `3*(A-1)+a`, column `3*(B-1)+b`
      real(dp), intent(in) :: h4(:, :, :, :)
      real(dp), intent(out) :: matrix(:, :)

      integer :: ia, ib, a, b

      do ib = 1, size(h4, 4)
         do ia = 1, size(h4, 3)
            do b = 1, 3
               do a = 1, 3
                  matrix(3*(ia - 1) + a, 3*(ib - 1) + b) = h4(a, b, ia, ib)
               end do
            end do
         end do
      end do
   end subroutine to_matrix

   subroutine differences_case(error, basis, aux_basis)
      !! The analytic Hessian against a five-point stencil of the fitted gradient
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: basis, aux_basis

      real(dp), parameter :: SHIFT(4) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]
      real(dp), parameter :: WEIGHT(4) = [-1.0_dp, 8.0_dp, -8.0_dp, 1.0_dp]
      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: hess4(:, :, :, :)
      real(dp) :: geo(3, N_ATOMS), g(N_COORD, 4), fd(N_COORD, N_COORD)
      real(dp) :: analytic(N_COORD, N_COORD)
      real(dp) :: dev
      integer :: y, k, atom, comp

      call converge(WATER, basis, aux_basis, mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge: "//err%get_full_trace())
      if (allocated(error)) return
      call df_rhf_hessian(mol, aux, scf%density, scf%orbitals, scf%orbital_energies, &
                          N_OCC, hess4, err, &
                          controls=response_controls_t(tol=1.0e-13_dp))
      call check(error,.not. err%has_error(), "the Hessian must run: "//err%get_full_trace())
      if (allocated(error)) return
      call to_matrix(hess4, analytic)

      call check(error, maxval(abs(analytic - transpose(analytic))) < 1.0e-9_dp, &
                 "the fitted Hessian is not symmetric")
      if (allocated(error)) return

      do y = 1, N_COORD
         atom = (y - 1)/3 + 1
         comp = y - 3*(atom - 1)
         do k = 1, 4
            geo = WATER
            geo(comp, atom) = geo(comp, atom) + SHIFT(k)*STEP
            call fitted_gradient(geo, basis, aux_basis, g(:, k), err)
            call check(error,.not. err%has_error(), "the displaced gradient failed")
            if (allocated(error)) return
         end do
         fd(:, y) = (WEIGHT(1)*g(:, 1) + WEIGHT(2)*g(:, 2) + WEIGHT(3)*g(:, 3) &
                     + WEIGHT(4)*g(:, 4))/(12.0_dp*STEP)
      end do
      fd = 0.5_dp*(fd + transpose(fd))
      dev = maxval(abs(analytic - fd))
      call check(error, dev < TOL_FD, "the fitted Hessian is not the derivative of the "// &
                 "fitted gradient")
      if (allocated(error)) write (*, "(a,es10.3)") "largest deviation ", dev
      call mol%destroy()
      call aux%destroy()
   end subroutine differences_case

   subroutine against_gradient_differences(error)
      !! s and p functions only, fitted in a set made for Coulomb and exchange
      type(error_type), allocatable, intent(out) :: error

      call differences_case(error, "6-31g", "cc-pvtz-jkfit")
   end subroutine against_gradient_differences

   subroutine against_gradient_differences_polarised(error)
      !! d functions on oxygen, and a different auxiliary family
      type(error_type), allocatable, intent(out) :: error

      call differences_case(error, "cc-pvdz", "cc-pvdz-rifit")
   end subroutine against_gradient_differences_polarised

   subroutine differs_from_exact(error)
      !! A fitted Hessian that equalled the exact one would be the exact one
      !! with fitted orbitals, which is the mistake of differentiating the
      !! wrong energy. The two are a few 1e-4 apart here.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: mol, aux
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: fitted(:, :, :, :), exact(:, :, :, :)

      call converge(WATER, "6-31g", "cc-pvdz-rifit", mol, aux, scf, err)
      call check(error,.not. err%has_error(), "the SCF must converge")
      if (allocated(error)) return
      call df_rhf_hessian(mol, aux, scf%density, scf%orbitals, scf%orbital_energies, &
                          N_OCC, fitted, err)
      call check(error,.not. err%has_error(), "the fitted Hessian must run")
      if (allocated(error)) return
      call rhf_hessian(mol, scf%density, scf%orbitals, scf%orbital_energies, N_OCC, &
                       exact, err)
      call check(error,.not. err%has_error(), "the exact Hessian must run")
      if (allocated(error)) return

      call check(error, maxval(abs(fitted - exact)) > MIN_FIT_EFFECT, &
                 "the fitted Hessian is indistinguishable from the exact one")
      call mol%destroy()
      call aux%destroy()
   end subroutine differs_from_exact

   subroutine metric_is_whole(error)
      !! The auxiliary sets used above invert without dropping a mode, which is
      !! the condition under which the analytic route is taken at all.
      type(error_type), allocatable, intent(out) :: error

      type(czt_molecule_t) :: aux
      type(error_t) :: err
      integer :: n_dropped

      call build_czt_molecule(NUMBERS, SYM, WATER, "cc-pvtz-jkfit", aux, err)
      call check(error,.not. err%has_error(), "the auxiliary basis must build")
      if (allocated(error)) return
      call df_metric_dropped(aux, n_dropped, err)
      call check(error,.not. err%has_error(), "the metric must factor")
      if (allocated(error)) return
      call check(error, n_dropped == 0, "cc-pvtz-jkfit on water should drop no mode")
      call aux%destroy()
   end subroutine metric_is_whole

end module test_mqc_czt_df_hessian

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_df_hessian, only: collect_mqc_czt_df_hessian_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_df_hessian", collect_mqc_czt_df_hessian_tests)]
   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
