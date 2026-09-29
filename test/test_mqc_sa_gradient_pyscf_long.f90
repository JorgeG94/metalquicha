!! SA-CASSCF per-root gradients against PySCF
module test_mqc_sa_gradient_pyscf_long
   !! Every root's analytic gradient of SA-2-CAS(2,2) against PySCF 2.14's
   !! `pyscf.grad.sacasscf`, fed this repository's basis JSON
   !! (`tools/sa_casscf/pyscf_ref.py --grad`; 6-31G* Cartesian, singlet CI).
   !! Measured agreement is 3e-10 to 5e-9 Hartree/Bohr; the 1e-7 tolerance
   !! leaves room for PySCF's convergence floor. The root energies are
   !! checked first, so a match cannot come from a different SA solution.
   !!
   !! Geometries: LiH as in `test_mqc_sa_gradient`, and
   !! `tools/sa_casscf/c2h4_planar.xyz` / `c2h4_twisted.xyz`, converted from
   !! Angstrom with `BOHR_TO_ANGSTROM` (PySCF uses the same value).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   use mqc_czt_sa_gradient, only: czt_sa_casscf_gradient
   implicit none
   private

   public :: collect_mqc_sa_gradient_pyscf_long_tests

   real(dp), parameter :: WEIGHTS(2) = [0.5_dp, 0.5_dp]
   real(dp), parameter :: GRAD_TOL = 1.0e-7_dp
   real(dp), parameter :: ENERGY_TOL = 1.0e-8_dp

   real(dp), parameter :: LIH_XYZ(3, 2) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                                   0.0_dp, 0.0_dp, 1.5949_dp], [3, 2])
   real(dp), parameter :: PLANAR_XYZ(3, 6) = reshape([ &
                                                     0.0_dp, 0.0_dp, 0.6695_dp, 0.0_dp, 0.0_dp, -0.6695_dp, &
                                                     0.0_dp, 0.9290_dp, 1.2320_dp, 0.0_dp, -0.9290_dp, 1.2320_dp, &
                                                   0.0_dp, 0.9290_dp, -1.2320_dp, 0.0_dp, -0.9290_dp, -1.2320_dp], [3, 6])
   real(dp), parameter :: TWISTED_XYZ(3, 6) = reshape([ &
                                                      0.0_dp, 0.0_dp, 0.6695_dp, 0.0_dp, 0.0_dp, -0.6695_dp, &
                                                      0.0_dp, 0.9290_dp, 1.2320_dp, 0.0_dp, -0.9290_dp, 1.2320_dp, &
                                                   0.9290_dp, 0.1500_dp, -1.2320_dp, -0.9290_dp, 0.1500_dp, -1.2320_dp], &
                                                      [3, 6])
   integer, parameter :: C2H4_Z(6) = [6, 6, 1, 1, 1, 1]
   character(len=2), parameter :: C2H4_SYM(6) = ["C ", "C ", "H ", "H ", "H ", "H "]

   real(dp), parameter :: LIH_E(2) = [-7.854944206511_dp, -7.725818481002_dp]
   real(dp), parameter :: LIH_G1(3, 2) = reshape([ &
                                                 -3.619309775677e-18_dp, 1.499672453459e-17_dp, -1.648036606632e-02_dp, &
                                            3.619309775677e-18_dp, -1.499672453459e-17_dp, 1.648036606632e-02_dp], [3, 2])
   real(dp), parameter :: LIH_G2(3, 2) = reshape([ &
                                                 -2.869140453094e-18_dp, 1.959597546003e-17_dp, 1.604264862932e-02_dp, &
                                           2.869140453094e-18_dp, -1.959597546003e-17_dp, -1.604264862932e-02_dp], [3, 2])
   real(dp), parameter :: PLANAR_E(2) = [-78.049709783650_dp, -77.673328496147_dp]
   real(dp), parameter :: PLANAR_G1(3, 6) = reshape([ &
                                                   8.224266709122e-16_dp, 6.224306469757e-15_dp, -9.094150559131e-03_dp, &
                                                    1.043767628118e-15_dp, 8.649217092947e-16_dp, 9.094150559129e-03_dp, &
                                                   -2.641710029602e-15_dp, 7.834521822513e-03_dp, 2.785171789867e-03_dp, &
                                                   1.758709370978e-15_dp, -7.834521822516e-03_dp, 2.785171789870e-03_dp, &
                                                   1.941736810178e-15_dp, 7.834521822542e-03_dp, -2.785171789871e-03_dp, &
                                          -2.924930450585e-15_dp, -7.834521822545e-03_dp, -2.785171789869e-03_dp], [3, 6])
   real(dp), parameter :: PLANAR_G2(3, 6) = reshape([ &
                                                  -1.893981551506e-15_dp, 4.376156540505e-15_dp, -1.398564672968e-01_dp, &
                                                   -2.167901380971e-15_dp, 2.153636819369e-15_dp, 1.398564672968e-01_dp, &
                                                    6.706692164267e-15_dp, 1.122750773920e-02_dp, 3.327531403611e-03_dp, &
                                                  -4.771019027572e-15_dp, -1.122750773920e-02_dp, 3.327531403613e-03_dp, &
                                                  -3.376804932814e-15_dp, 1.122750773917e-02_dp, -3.327531403595e-03_dp, &
                                           5.503014728595e-15_dp, -1.122750773917e-02_dp, -3.327531403594e-03_dp], [3, 6])
   real(dp), parameter :: TWISTED_E(2) = [-77.920166075657_dp, -77.807933612709_dp]
   real(dp), parameter :: TWISTED_G1(3, 6) = reshape([ &
                                                  -1.336915876781e-15_dp, 4.102368389150e-03_dp, -1.073977872407e-01_dp, &
                                                   8.768678441281e-16_dp, -1.063783052716e-02_dp, 1.148954013387e-01_dp, &
                                                  -1.363227042769e-16_dp, 6.587620194191e-03_dp, -1.751043735271e-03_dp, &
                                                 -2.740066268959e-16_dp, -7.584541215597e-03_dp, -1.713777059795e-03_dp, &
                                                   1.104349887899e-02_dp, 3.766191579705e-03_dp, -2.016396651427e-03_dp, &
                                           -1.104349887899e-02_dp, 3.766191579704e-03_dp, -2.016396651425e-03_dp], [3, 6])
   real(dp), parameter :: TWISTED_G2(3, 6) = reshape([ &
                                                   1.748567071356e-15_dp, -5.338834546566e-03_dp, 7.842853882863e-03_dp, &
                                                   -4.627770003525e-15_dp, 1.201630622033e-02_dp, 5.509225079190e-02_dp, &
                                                  1.001058894900e-16_dp, -2.412719525761e-04_dp, -2.576978154404e-02_dp, &
                                                  9.709418650469e-17_dp, -8.386996045609e-04_dp, -2.174351633600e-02_dp, &
                                                  8.612767917740e-03_dp, -2.798750058319e-03_dp, -7.710903397360e-03_dp, &
                                          -8.612767917737e-03_dp, -2.798750058318e-03_dp, -7.710903397360e-03_dp], [3, 6])

contains

   subroutine collect_mqc_sa_gradient_pyscf_long_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lih_sto3g_sa2_matches_pyscf", test_lih), &
                  new_unittest("c2h4_planar_sa2_matches_pyscf", test_planar), &
                  new_unittest("c2h4_twisted_sa2_matches_pyscf", test_twisted) &
                  ]
   end subroutine collect_mqc_sa_gradient_pyscf_long_tests

   subroutine compare(z, sym, xyz_angstrom, basis, n_electrons, n_inactive, e_ref, &
                      g1_ref, g2_ref, error)
      !! SA-2-CAS(2,2) at `xyz_angstrom`, both roots' energies and gradients
      !! against the PySCF numbers
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz_angstrom(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_electrons, n_inactive
      real(dp), intent(in) :: e_ref(2), g1_ref(:, :), g2_ref(:, :)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: result
      real(dp), allocatable :: g(:, :)
      real(dp) :: worst
      integer :: root

      call build_czt_molecule(z, sym, xyz_angstrom/BOHR_TO_ANGSTROM, basis, mol, err)
      call run_czt_rhf(mol, n_electrons, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "RHF should converge")
      if (allocated(error)) return
      call run_czt_casscf(mol, scf%orbitals, n_inactive, 2, 1, 1, result, err, &
                          max_iterations=400, gradient_tol=1.0e-10_dp, n_states=2, &
                          weights=WEIGHTS)
      call check(error,.not. err%has_error() .and. result%converged, &
                 "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return
      call check(error, maxval(abs(result%energies - e_ref)) < ENERGY_TOL, &
                 "both root energies should match PySCF's")
      if (allocated(error)) return

      do root = 1, 2
         call czt_sa_casscf_gradient(mol, result%orbitals, n_inactive, 2, 1, 1, &
                                     result%ci_vectors, result%energies, WEIGHTS, root, &
                                     g, err)
         call check(error,.not. err%has_error(), "the root gradient should build")
         if (allocated(error)) return
         if (root == 1) then
            worst = maxval(abs(g - g1_ref))
         else
            worst = maxval(abs(g - g2_ref))
         end if
         write (*, "(a,i0,a,es10.2)") "    root ", root, " max |mqc - PySCF| = ", worst
         call check(error, worst < GRAD_TOL, "the root gradient should match PySCF's")
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine compare

   subroutine test_lih(error)
      !! LiH/STO-3G
      type(error_type), allocatable, intent(out) :: error
      call compare([3, 1], ["Li", "H "], LIH_XYZ, "sto-3g", 4, 1, LIH_E, LIH_G1, LIH_G2, error)
   end subroutine test_lih

   subroutine test_planar(error)
      !! Planar C2H4/6-31G*: S0 and the ionic V state
      type(error_type), allocatable, intent(out) :: error
      call compare(C2H4_Z, C2H4_SYM, PLANAR_XYZ, "6-31g*", 16, 7, PLANAR_E, PLANAR_G1, &
                   PLANAR_G2, error)
   end subroutine test_planar

   subroutine test_twisted(error)
      !! Twisted, pyramidalised C2H4/6-31G*: S0 and S1 3 eV apart
      type(error_type), allocatable, intent(out) :: error
      call compare(C2H4_Z, C2H4_SYM, TWISTED_XYZ, "6-31g*", 16, 7, TWISTED_E, TWISTED_G1, &
                   TWISTED_G2, error)
   end subroutine test_twisted

end module test_mqc_sa_gradient_pyscf_long

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_gradient_pyscf_long, only: collect_mqc_sa_gradient_pyscf_long_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_gradient_pyscf_long", &
                               collect_mqc_sa_gradient_pyscf_long_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
