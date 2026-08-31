!! The double hybrid's dipole derivative, against a difference of its dipole
module test_mqc_dh_dipole_deriv
   !! What a double-hybrid infrared intensity needs beyond the reference's
   !!
   !! A double hybrid's dipole is the *relaxed* one: the perturbative term
   !! contributes through the MP2 relaxed density, and on water/6-31G that
   !! contribution is 1.26e-02 -- small against a dipole of order one, and far
   !! too large to leave out of an intensity.
   !!
   !! **This settles a reading, not just a number.** `mp2_perturbed_response`
   !! produces `ddrel`, and `correlation_dipole_derivatives` takes it as the AO
   !! derivative expressed in the MO frame -- `dD_AO/dR = C ddrel C^T`, with no
   !! orbital-rotation terms added. Nothing in that routine's own code says so:
   !! `ddrel` is contracted there against a skeleton first-order Fock, which
   !! would be equally consistent with a convention that kept the rotations
   !! separate. Getting it wrong leaves a derivative that is smooth, has the
   !! right sign and obeys the sum rule, and is wrong by the size of the
   !! response. Differencing the dipole is what tells the two apart.
   !!
   !! The dipole being differenced is itself pinned: `-c Tr(C dm1mo C^T r)`
   !! agrees with a field-differenced PySCF B2PLYP energy to 1.8e-08, which is
   !! that difference's own precision. So this rests on a quantity with an
   !! external reference rather than on its own internal consistency.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_libcint_integrals, only: libcint_molecule_t, build_libcint_molecule
   use mqc_libcint_rhf, only: rhf_result_t, run_libcint_rhf
   use mqc_libcint_mp2_gradient, only: libcint_mp2_gradient
   use mqc_libcint_mp2_hessian, only: mp2_correlation_hessian
   use mqc_libcint_multipole, only: multipole_matrices
   use mqc_libcint_xc, only: xc_context_t, xc_context_create, xc_available
   use omp_lib, only: omp_set_num_threads, omp_get_max_threads
   implicit none
   private

   public :: collect_mqc_dh_dipole_deriv_tests

   integer, parameter :: WATER_Z(3) = [8, 1, 1]
   character(len=2), parameter :: WATER_SYM(3) = ["O ", "H ", "H "]
   real(dp), parameter :: WATER(3, 3) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 1.814137_dp, &
                           0.0_dp, 1.756_dp, -0.4543_dp], [3, 3])
   integer, parameter :: WATER_NELEC = 10
   character(len=*), parameter :: BASIS = "6-31g"
   character(len=*), parameter :: FUNCTIONAL = "b2plyp"
   integer, parameter :: GRID_LEVEL = 5

   !> One column, as the rest of this ladder takes one.
   integer, parameter :: PERT_ATOM = 1
   integer, parameter :: PERT_CART = 3
   real(dp), parameter :: STEP = 2.0e-3_dp

   !> Fifteen times the measured worst, 6.687e-10 -- the stencil's own floor.
   !>
   !> Worth stating what this number was before the orbital rotations went in:
   !> **1.110e-02**, against a correlation dipole of 1.26e-02. Seven orders,
   !> from two terms whose absence left a derivative that was smooth, correctly
   !> signed, and summed to the molecular charge.
   real(dp), parameter :: TOL = 1.0e-8_dp

   !> `-c Tr(C dm1mo C^T r)` at the reference geometry, about the nuclear
   !> centroid, from a field difference of PySCF's own B2PLYP energy:
   !> `E(F) = E0 - mu.F` with the field in `hcore`, so this is the electronic
   !> part and the nuclei cancel out of the Kohn-Sham/double-hybrid difference.
   real(dp), parameter :: PYSCF_MU_CORR(3) = [0.0_dp, -1.2571831434e-02_dp, &
                                              -9.7393535725e-03_dp]

contains

   subroutine collect_mqc_dh_dipole_deriv_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("correlation_dipole_against_a_field_difference", &
                               test_dipole), &
                  new_unittest("correlation_derivative_differences_its_dipole", &
                               test_derivative), &
                  new_unittest("hartree_fock_reference_is_refused", test_refusal) &
                  ]
   end subroutine collect_mqc_dh_dipole_deriv_tests

   !> The nuclear centroid, which both halves of the assembly expand about
   pure function centroid(coords) result(o)
      real(dp), intent(in) :: coords(3, 3)
      real(dp) :: o(3)

      o = (coords(:, 1) + coords(:, 2) + coords(:, 3))/3.0_dp
   end function centroid

   !> `c mu_corr` at one geometry, about a fixed origin
   subroutine mu_corr_at(coords, origin, mu, err, ok)
      real(dp), intent(in) :: coords(3, 3), origin(3)
      real(dp), intent(out) :: mu(3)
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      real(dp), allocatable :: grad(:, :), dm1mo(:, :), dip(:, :, :), drel(:, :)
      integer :: a

      ok = .false.
      call build_libcint_molecule(WATER_Z, WATER_SYM, coords, BASIS, mol, err)
      if (err%has_error()) return
      call xc_context_create(mol, FUNCTIONAL, ctx, err, level=GRID_LEVEL)
      if (.not. err%has_error()) &
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                              .false., scf, err, xc=ctx)
      if (.not. err%has_error()) &
         call libcint_mp2_gradient(mol, scf%orbitals, scf%orbital_energies, &
                                   WATER_NELEC/2, grad, err, xc=ctx, &
                                   scf_density=scf%density, &
                                   pt2_scale=ctx%pt2_fraction, &
                                   relaxed_density_mo=dm1mo)
      if (.not. err%has_error()) call multipole_matrices(mol, origin, 1, dip, err)
      if (err%has_error()) then
         call mol%destroy()
         return
      end if

      drel = matmul(scf%orbitals, matmul(dm1mo, transpose(scf%orbitals)))
      do a = 1, 3
         mu(a) = -ctx%pt2_fraction*sum(drel*dip(:, :, a))
      end do
      call mol%destroy()
      ok = .true.
   end subroutine mu_corr_at

   subroutine test_dipole(error)
      !! The quantity the derivative below is differenced against
      type(error_type), allocatable, intent(out) :: error

      real(dp) :: mu(3)
      type(error_t) :: err
      logical :: ok

      if (.not. xc_available()) return
      call mu_corr_at(WATER, centroid(WATER), mu, err, ok)
      call check(error, ok, "the correlation dipole did not evaluate: "// &
                 err%get_message())
      if (allocated(error)) return

      write (*, "(a, es10.3)") "        max |ours - field difference| = ", &
         maxval(abs(mu - PYSCF_MU_CORR))
      ! Ten times the measured 1.79e-08, which is the field difference's own
      ! floor rather than ours: PySCF's energy there is converged to 1e-14 and
      ! divided by a step of 1e-4.
      call check(error, maxval(abs(mu - PYSCF_MU_CORR)) < 2.0e-7_dp, &
                 "the perturbative term's dipole disagrees with a field "// &
                 "difference of PySCF's own double-hybrid energy")
   end subroutine test_dipole

   subroutine test_derivative(error)
      !! One column, against the 7-point stencil
      type(error_type), allocatable, intent(out) :: error

      real(dp), parameter :: WEIGHT(6) = [-1.0_dp, 9.0_dp, -45.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]
      integer, parameter :: OFFSET(6) = [-3, -2, -1, 1, 2, 3]
      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(xc_context_t) :: ctx
      type(error_t) :: err
      real(dp), allocatable :: hc(:, :, :, :), hr(:, :, :, :), ddip(:, :)
      real(dp) :: origin(3), coords(3, 3), mu(3), fd(3)
      logical :: ok
      integer :: point, threads, col

      if (.not. xc_available()) return
      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      origin = centroid(WATER)

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, BASIS, mol, err)
      if (.not. err%has_error()) &
         call xc_context_create(mol, FUNCTIONAL, ctx, err, level=GRID_LEVEL)
      if (.not. err%has_error()) &
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-14_dp, 1.0e-12_dp, &
                              .false., scf, err, xc=ctx)
      if (.not. err%has_error()) &
         call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%density, WATER_NELEC/2, 0, hc, hr, err, &
                                      xc=ctx, pt2_scale=ctx%pt2_fraction, &
                                      dipole_derivatives=ddip)
      call mol%destroy()
      call check(error,.not. err%has_error(), &
                 "the correlation dipole derivative did not evaluate: "// &
                 err%get_message())
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if

      fd = 0.0_dp
      do point = 1, 6
         coords = WATER
         coords(PERT_CART, PERT_ATOM) = coords(PERT_CART, PERT_ATOM) &
                                        + real(OFFSET(point), dp)*STEP
         call mu_corr_at(coords, origin, mu, err, ok)
         if (.not. ok) exit
         fd = fd + WEIGHT(point)*mu
      end do
      call omp_set_num_threads(threads)
      call check(error, ok, "a displaced correlation dipole failed")
      if (allocated(error)) return
      fd = fd/(60.0_dp*STEP)

      col = 3*(PERT_ATOM - 1) + PERT_CART
      write (*, "(a, es10.3)") "        max |analytic - fd| = ", &
         maxval(abs(ddip(:, col) - fd))
      call check(error, maxval(abs(ddip(:, col) - fd)) < TOL, &
                 "the perturbative term's dipole derivative does not difference "// &
                 "its own dipole")
   end subroutine test_derivative

   subroutine test_refusal(error)
      !! Asking for it over a Hartree-Fock reference is refused, not approximated
      type(error_type), allocatable, intent(out) :: error

      type(libcint_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(error_t) :: err
      real(dp), allocatable :: hc(:, :, :, :), hr(:, :, :, :), ddip(:, :)

      call build_libcint_molecule(WATER_Z, WATER_SYM, WATER, "sto-3g", mol, err)
      if (.not. err%has_error()) &
         call run_libcint_rhf(mol, WATER_NELEC, 300, 1.0e-12_dp, 1.0e-10_dp, &
                              .false., scf, err)
      if (.not. err%has_error()) &
         call mp2_correlation_hessian(mol, scf%orbitals, scf%orbital_energies, &
                                      scf%density, WATER_NELEC/2, 0, hc, hr, err, &
                                      dipole_derivatives=ddip)
      call mol%destroy()
      call check(error, err%has_error(), &
                 "a dipole derivative over a Hartree-Fock reference should be "// &
                 "refused rather than returned")
   end subroutine test_refusal

end module test_mqc_dh_dipole_deriv

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_dh_dipole_deriv, only: collect_mqc_dh_dipole_deriv_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_dh_dipole_deriv", collect_mqc_dh_dipole_deriv_tests)]
   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do
   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
