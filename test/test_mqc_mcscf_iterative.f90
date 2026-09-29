!! The iterative CASSCF Newton step against the explicit one
module test_mqc_mcscf_iterative
   !! `run_czt_casscf` with `iterative_hessian` forced on and off must reach
   !! the same energy in the same number of macro-iterations: the iterative
   !! step (`iterative_newton_step`, a Krylov subspace over Hessian-vector
   !! products) and the explicit one (`orbital_hessian` diagonalised) share
   !! `level_shifted_step`, so they differ only in how well the subspace
   !! resolves the Newton equations. Water/STO-3G CAS(6,5) is the
   !! symmetry-trapped case whose saddle escape the explicit path is known to
   !! need; H2O/6-31G CAS(4,4) SA-3 covers state averaging.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   implicit none
   private

   public :: collect_mqc_mcscf_iterative_tests

   real(dp), parameter :: WATER_EQ(3, 3) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                                    0.0_dp, -0.7572_dp, 0.5865_dp, &
                                                    0.0_dp, 0.7572_dp, 0.5865_dp], [3, 3])
   real(dp), parameter :: WATER_ASYM(3, 3) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                                      0.0_dp, -0.7572_dp, 0.5865_dp, &
                                                      0.0_dp, 0.8100_dp, 0.6200_dp], [3, 3])

contains

   subroutine collect_mqc_mcscf_iterative_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("water_sto3g_cas65_iterative_matches_explicit", test_water), &
                  new_unittest("h2o_631g_cas44_sa3_iterative_matches_explicit", test_h2o_sa3) &
                  ]
   end subroutine collect_mqc_mcscf_iterative_tests

   subroutine compare(xyz, basis, n_inactive, n_active, n_half, weights, error)
      !! Converge both ways from the same RHF orbitals and compare
      real(dp), intent(in) :: xyz(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_inactive, n_active, n_half
      real(dp), intent(in) :: weights(:)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: explicit, iterative
      integer :: n_states

      n_states = size(weights)
      call build_czt_molecule([8, 1, 1], ["O ", "H ", "H "], xyz/BOHR_TO_ANGSTROM, basis, mol, err)
      call run_czt_rhf(mol, 10, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "RHF should converge")
      if (allocated(error)) return

      if (n_states > 1) then
         call run_czt_casscf(mol, scf%orbitals, n_inactive, n_active, n_half, n_half, explicit, &
                             err, max_iterations=200, gradient_tol=1.0e-8_dp, n_states=n_states, &
                             weights=weights, iterative_hessian=.false.)
         call run_czt_casscf(mol, scf%orbitals, n_inactive, n_active, n_half, n_half, iterative, &
                             err, max_iterations=200, gradient_tol=1.0e-8_dp, n_states=n_states, &
                             weights=weights, iterative_hessian=.true.)
      else
         call run_czt_casscf(mol, scf%orbitals, n_inactive, n_active, n_half, n_half, explicit, &
                             err, max_iterations=200, gradient_tol=1.0e-8_dp, &
                             iterative_hessian=.false.)
         call run_czt_casscf(mol, scf%orbitals, n_inactive, n_active, n_half, n_half, iterative, &
                             err, max_iterations=200, gradient_tol=1.0e-8_dp, &
                             iterative_hessian=.true.)
      end if
      call check(error,.not. err%has_error(), "both CASSCF runs should succeed")
      if (allocated(error)) return
      call check(error, explicit%converged .and. iterative%converged, "both should converge")
      if (allocated(error)) return
      call check(error, iterative%energy, explicit%energy, &
                 "the iterative step should reach the explicit step's energy", thr=1.0e-10_dp)
      if (allocated(error)) return
      call check(error, iterative%iterations, explicit%iterations, &
                 "and in the same number of macro-iterations")
      call mol%destroy()
   end subroutine compare

   subroutine test_water(error)
      !! Water/STO-3G CAS(6,5), from symmetric RHF orbitals
      type(error_type), allocatable, intent(out) :: error
      call compare(WATER_EQ, "sto-3g", 2, 5, 3, [1.0_dp], error)
   end subroutine test_water

   subroutine test_h2o_sa3(error)
      !! H2O (one bond stretched)/6-31G CAS(4,4) SA-3
      type(error_type), allocatable, intent(out) :: error
      call compare(WATER_ASYM, "6-31g", 3, 4, 2, [1.0_dp, 1.0_dp, 1.0_dp]/3.0_dp, error)
   end subroutine test_h2o_sa3

end module test_mqc_mcscf_iterative

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_mcscf_iterative, only: collect_mqc_mcscf_iterative_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_mcscf_iterative", collect_mqc_mcscf_iterative_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
