!! SA-CASSCF nonadiabatic couplings against PySCF
module test_mqc_sa_nac_pyscf_long
   !! `czt_sa_casscf_nacs` against PySCF 2.14's
   !! `pyscf.nac.sacasscf.NonAdiabaticCouplings`, fed this repository's basis
   !! JSON (`tools/sa_casscf/pyscf_ref.py --nac`): twisted and planar C2H4
   !! 6-31G* SA-2-CAS(2,2), and H2O/6-31G CAS(4,4) SA-3 for all three pairs.
   !! Not linear LiH: its active space holds one member of a degenerate pi
   !! pair, so the coupling has a null orbital direction and PySCF lands
   !! slightly off the symmetric solution (x/y components of 1e-4 that
   !! symmetry makes zero). And H2O with one O-H bond stretched, since at the
   !! symmetric geometry SA-3 breaks C2v into two mirror solutions and a run
   !! can land on either.
   !!
   !! `d_IJ` with the CSF term is PySCF's `use_etfs=False`, without it
   !! `use_etfs=True`. `h_IJ = (E_J - E_I) d_IJ` is the negative of PySCF's
   !! `mult_ediff=True`, which scales by `E_I - E_J`. The CI phase makes each
   !! pair's sign arbitrary: one sign is chosen per pair, from `d_IJ` with
   !! the CSF term, and applied to every quantity of that pair.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   use mqc_czt_sa_nac, only: czt_sa_casscf_nacs
   implicit none
   private

   public :: collect_mqc_sa_nac_pyscf_long_tests

   real(dp), parameter :: NAC_TOL = 2.0e-6_dp
      !! Measured: C2H4 1e-8 to 2e-7, H2O 2e-7 to 9e-7 (pair (2,3), whose
      !! |d| is about 2). PySCF's SA-CASSCF orbitals stop near 1e-7, and a
      !! NAC divides by an energy gap.

   real(dp), parameter :: H2O_XYZ(3, 3) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                                   0.0_dp, -0.7572_dp, 0.5865_dp, &
                                                   0.0_dp, 0.8100_dp, 0.6200_dp], [3, 3])
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

   real(dp), parameter :: PLANAR_D_CSF(3, 6) = reshape([ &
                                            -2.91939125846388e-15_dp, -8.83084132032506e-16_dp, 2.93691255700986e-01_dp, &
                                             2.82547819160195e-15_dp, -1.02663141272273e-14_dp, 2.93691255701334e-01_dp, &
                                            1.10992115869798e-15_dp, -3.27514055158010e-03_dp, -1.17929694881196e-02_dp, &
                                            -3.34094916786072e-15_dp, 3.27514055158273e-03_dp, -1.17929694881110e-02_dp, &
                                             6.07386742854449e-15_dp, 3.27514055156824e-03_dp, -1.17929694881150e-02_dp, &
                                    -3.95464475147339e-15_dp, -3.27514055156170e-03_dp, -1.17929694881146e-02_dp], [3, 6])
   real(dp), parameter :: PLANAR_D_NOCSF(3, 6) = reshape([ &
                                             -3.20152135629434e-15_dp, 2.44705432940720e-14_dp, 2.35859389762693e-02_dp, &
                                             2.86960123198889e-15_dp, -2.25220243040815e-14_dp, 2.35859389761406e-02_dp, &
                                            8.97218988805637e-16_dp, -3.27514055157597e-03_dp, -1.17929694881171e-02_dp, &
                                            -3.04919385240758e-15_dp, 3.27514055157007e-03_dp, -1.17929694881095e-02_dp, &
                                             5.92172639946472e-15_dp, 3.27514055156800e-03_dp, -1.17929694881007e-02_dp, &
                                    -3.43783141155704e-15_dp, -3.27514055156531e-03_dp, -1.17929694881032e-02_dp], [3, 6])
   real(dp), parameter :: PLANAR_H_CSF(3, 6) = reshape([ &
                                            -1.09880424056408e-15_dp, -2.65263756704995e-15_dp, 1.10539892949120e-01_dp, &
                                             1.06345711968565e-15_dp, -2.63649917048368e-15_dp, 1.10539892949266e-01_dp, &
                                            4.17753563120655e-16_dp, -1.23270161755860e-03_dp, -4.43865303942421e-03_dp, &
                                            -1.25747075770455e-15_dp, 1.23270161756144e-03_dp, -4.43865303941926e-03_dp, &
                                             2.28609003448543e-15_dp, 1.23270161755309e-03_dp, -4.43865303941992e-03_dp, &
                                    -1.48845427488480e-15_dp, -1.23270161755145e-03_dp, -4.43865303941820e-03_dp], [3, 6])
   real(dp), parameter :: TWISTED_D_CSF(3, 6) = reshape([ &
                                            -1.00368224472206e-02_dp, -2.86704446635357e-14_dp, 4.36309782629483e-13_dp, &
                                            -2.89759255709431e-02_dp, 6.42731139249393e-14_dp, -4.50496366395984e-13_dp, &
                                            -3.00407713389077e-01_dp, -2.56206894622395e-14_dp, 1.29531957810060e-15_dp, &
                                            2.88180234439744e-01_dp, -9.32954689908982e-15_dp, -2.66870270236329e-14_dp, &
                                            2.41248877186812e-02_dp, -4.09079089534023e-01_dp, -4.62005684119439e-02_dp, &
                                       2.41248877186191e-02_dp, 4.09079089534031e-01_dp, 4.62005684120703e-02_dp], [3, 6])
   real(dp), parameter :: TWISTED_D_NOCSF(3, 6) = reshape([ &
                                            -9.25470428431194e-03_dp, 1.34029740169407e-14_dp, -1.81506108167227e-13_dp, &
                                             -2.16481331713097e-02_dp, 1.73328676387703e-14_dp, 2.39206618746531e-13_dp, &
                                           -3.58530280102885e-01_dp, -7.51899872461915e-14_dp, -1.45367660136363e-14_dp, &
                                             3.43231133726456e-01_dp, 4.39319308806693e-14_dp, -4.08376203476878e-14_dp, &
                                            2.31009919160520e-02_dp, -3.57411718637522e-01_dp, -4.05865702584664e-02_dp, &
                                       2.31009919160096e-02_dp, 3.57411718637528e-01_dp, 4.05865702585275e-02_dp], [3, 6])
   real(dp), parameter :: TWISTED_H_CSF(3, 6) = reshape([ &
                                             -1.12645730340505e-03_dp, 4.62444476570641e-15_dp, 1.21605011090426e-14_dp, &
                                           -3.25203949303129e-03_dp, -1.89070706235382e-15_dp, -8.42943992184971e-15_dp, &
                                             -3.37154975620955e-02_dp, 5.88268177850195e-16_dp, 3.25575028678269e-15_dp, &
                                            3.23431774839730e-02_dp, -3.95895331205478e-15_dp, -1.06775144057187e-15_dp, &
                                            2.70759556700434e-03_dp, -4.59119537587415e-02_dp, -5.18520358245817e-03_dp, &
                                       2.70759556699497e-03_dp, 4.59119537587434e-02_dp, 5.18520358246024e-03_dp], [3, 6])
   real(dp), parameter :: H2O_12_D_CSF(3, 3) = reshape([ &
                                             -1.30053028158446e-01_dp, 9.81080971567390e-17_dp, 4.87103323221034e-15_dp, &
                                            6.34220365886469e-02_dp, -2.25144827114702e-15_dp, -9.23101491280113e-16_dp, &
                                     1.81453868478506e-01_dp, -6.20386613238417e-15_dp, -2.77352182576870e-15_dp], [3, 3])
   real(dp), parameter :: H2O_12_D_NOCSF(3, 3) = reshape([ &
                                             -5.01967795993000e-01_dp, 2.91526447734173e-16_dp, 2.46369751010304e-16_dp, &
                                             1.75108297066698e-01_dp, -6.46301239730104e-15_dp, 4.93620524971522e-15_dp, &
                                     3.26859498926297e-01_dp, -1.64711903326469e-15_dp, -1.89053079742817e-15_dp], [3, 3])
   real(dp), parameter :: H2O_12_H_CSF(3, 3) = reshape([ &
                                            -3.50877102170385e-02_dp, 2.64691144252716e-17_dp, -2.46057586961616e-15_dp, &
                                             1.71109744440991e-02_dp, -2.60583185656870e-15_dp, 2.41548671263785e-15_dp, &
                                       4.89554210700700e-02_dp, 3.24626829767588e-16_dp, 1.02807337779702e-15_dp], [3, 3])
   real(dp), parameter :: H2O_13_D_CSF(3, 3) = reshape([ &
                                            -6.44378688934969e-02_dp, -5.46764426299973e-15_dp, 4.86057355307698e-15_dp, &
                                            -1.17403219035847e-01_dp, 3.35397998280212e-15_dp, -4.95319277490131e-15_dp, &
                                      1.48659815917096e-01_dp, 1.62710126052615e-15_dp, -2.65829576482406e-15_dp], [3, 3])
   real(dp), parameter :: H2O_13_D_NOCSF(3, 3) = reshape([ &
                                             3.59137275064012e-02_dp, -1.19253190920358e-15_dp, 4.71200000688934e-15_dp, &
                                            -3.17835019896583e-01_dp, 7.14394106701479e-15_dp, -6.10481102989695e-15_dp, &
                                     2.81921292390187e-01_dp, -9.51285933928503e-16_dp, -4.85734300684577e-15_dp], [3, 3])
   real(dp), parameter :: H2O_13_H_CSF(3, 3) = reshape([ &
                                            -2.28923656478299e-02_dp, -7.21204258074687e-16_dp, 3.94512429225099e-16_dp, &
                                            -4.17089743120360e-02_dp, 2.07972211070534e-15_dp, -1.75968420548961e-15_dp, &
                                      5.28132745783114e-02_dp, 3.56003639711750e-16_dp, -1.05541540070426e-15_dp], [3, 3])
   real(dp), parameter :: H2O_23_D_CSF(3, 3) = reshape([ &
                                             -1.21547199698510e-15_dp, 1.90376987290084e+00_dp, 9.54156713237280e-02_dp, &
                                             1.19766264886709e-15_dp, -9.66447686046437e-01_dp, 6.64882946162311e-01_dp, &
                                     3.59031776986858e-17_dp, -1.08755490323977e+00_dp, -7.55760764755975e-01_dp], [3, 3])
   real(dp), parameter :: H2O_23_D_NOCSF(3, 3) = reshape([ &
                                             -1.18155155988382e-15_dp, 1.95488007774325e+00_dp, 9.79191606536840e-02_dp, &
                                             1.26491734568741e-15_dp, -9.16597817165234e-01_dp, 6.72113352519697e-01_dp, &
                                    -8.33657858059114e-17_dp, -1.03828226057803e+00_dp, -7.70032513173364e-01_dp], [3, 3])
   real(dp), parameter :: H2O_23_H_CSF(3, 3) = reshape([ &
                                             -1.03882989802487e-16_dp, 1.62709882873805e-01_dp, 8.15491038407347e-03_dp, &
                                             1.02360874309126e-16_dp, -8.25995788874777e-02_dp, 5.68256845718543e-02_dp, &
                                     3.06854446768969e-18_dp, -9.29502737929871e-02_dp, -6.45927573833736e-02_dp], [3, 3])

contains

   subroutine collect_mqc_sa_nac_pyscf_long_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("c2h4_planar_sa2_nac_matches_pyscf", test_planar), &
                  new_unittest("c2h4_twisted_sa2_nac_matches_pyscf", test_twisted), &
                  new_unittest("h2o_631g_cas44_sa3_nacs_match_pyscf", test_cas44) &
                  ]
   end subroutine collect_mqc_sa_nac_pyscf_long_tests

   subroutine converge(z, sym, xyz_angstrom, basis, n_electrons, n_inactive, n_active, &
                       n_half, weights, mol, result, error)
      !! A tightly converged SA-CASSCF at `xyz_angstrom`
      integer, intent(in) :: z(:)
      character(len=2), intent(in) :: sym(:)
      real(dp), intent(in) :: xyz_angstrom(:, :)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_electrons, n_inactive, n_active, n_half
      real(dp), intent(in) :: weights(:)
      type(czt_molecule_t), intent(out) :: mol
      type(casscf_result_t), intent(out) :: result
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(rhf_result_t) :: scf

      call build_czt_molecule(z, sym, xyz_angstrom/BOHR_TO_ANGSTROM, basis, mol, err)
      call run_czt_rhf(mol, n_electrons, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "RHF should converge")
      if (allocated(error)) return
      call run_czt_casscf(mol, scf%orbitals, n_inactive, n_active, n_half, n_half, result, &
                          err, max_iterations=500, gradient_tol=1.0e-10_dp, &
                          n_states=size(weights), weights=weights)
      call check(error,.not. err%has_error() .and. result%converged, &
                 "the SA-CASSCF should converge")
   end subroutine converge

   subroutine compare_pairs(mol, result, n_inactive, n_active, n_half, weights, pairs, &
                            d_csf, d_nocsf, h_csf, error)
      !! Every pair's `d_IJ` (with and without the CSF term) and `h_IJ` against
      !! PySCF, up to one sign per pair
      type(czt_molecule_t), intent(in) :: mol
      type(casscf_result_t), intent(in) :: result
      integer, intent(in) :: n_inactive, n_active, n_half
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: pairs(:, :)
      real(dp), intent(in) :: d_csf(:, :, :), d_nocsf(:, :, :), h_csf(:, :, :)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      real(dp), allocatable :: d(:, :, :), h(:, :, :), csf(:, :, :), ediff(:)
      real(dp), allocatable :: d0(:, :, :), h0(:, :, :)
      real(dp) :: sgn, worst_d, worst_h, worst_d0
      integer :: ip

      call czt_sa_casscf_nacs(mol, result%orbitals, n_inactive, n_active, n_half, n_half, &
                              result%ci_vectors, result%energies, weights, pairs, d, h, &
                              csf, ediff, err, include_csf=.true.)
      call check(error,.not. err%has_error(), "the NACs with the CSF term should build")
      if (allocated(error)) return
      call czt_sa_casscf_nacs(mol, result%orbitals, n_inactive, n_active, n_half, n_half, &
                              result%ci_vectors, result%energies, weights, pairs, d0, h0, &
                              csf, ediff, err, include_csf=.false.)
      call check(error,.not. err%has_error(), "the NACs without the CSF term should build")
      if (allocated(error)) return

      do ip = 1, size(pairs, 2)
         sgn = merge(1.0_dp, -1.0_dp, sum(d(:, :, ip)*d_csf(:, :, ip)) >= 0.0_dp)
         worst_d = maxval(abs(sgn*d(:, :, ip) - d_csf(:, :, ip)))
         worst_h = maxval(abs(sgn*h(:, :, ip) - h_csf(:, :, ip)))
         worst_d0 = maxval(abs(sgn*d0(:, :, ip) - d_nocsf(:, :, ip)))
         write (*, "(a,i0,a,i0,a,3es10.2)") "    pair (", pairs(1, ip), ",", pairs(2, ip), &
            ") max |mqc - PySCF|, d / h / d without CSF: ", worst_d, worst_h, worst_d0
         call check(error, max(worst_d, worst_h, worst_d0) < NAC_TOL, &
                    "the NAC should match PySCF's")
         if (allocated(error)) return
      end do
   end subroutine compare_pairs

   subroutine test_planar(error)
      !! Planar C2H4: S0 and the ionic V state
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(casscf_result_t) :: result
      real(dp), parameter :: W(2) = [0.5_dp, 0.5_dp]

      call converge(C2H4_Z, C2H4_SYM, PLANAR_XYZ, "6-31g*", 16, 7, 2, 1, W, mol, result, error)
      if (allocated(error)) return
      call compare_pairs(mol, result, 7, 2, 1, W, reshape([1, 2], [2, 1]), &
                         reshape(PLANAR_D_CSF, [3, 6, 1]), reshape(PLANAR_D_NOCSF, [3, 6, 1]), &
                         reshape(PLANAR_H_CSF, [3, 6, 1]), error)
      call mol%destroy()
   end subroutine test_planar

   subroutine test_twisted(error)
      !! Twisted, pyramidalised C2H4: S0 and S1 3 eV apart
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(casscf_result_t) :: result
      real(dp), parameter :: W(2) = [0.5_dp, 0.5_dp]

      call converge(C2H4_Z, C2H4_SYM, TWISTED_XYZ, "6-31g*", 16, 7, 2, 1, W, mol, result, error)
      if (allocated(error)) return
      call compare_pairs(mol, result, 7, 2, 1, W, reshape([1, 2], [2, 1]), &
                         reshape(TWISTED_D_CSF, [3, 6, 1]), reshape(TWISTED_D_NOCSF, [3, 6, 1]), &
                         reshape(TWISTED_H_CSF, [3, 6, 1]), error)
      call mol%destroy()
   end subroutine test_twisted

   subroutine test_cas44(error)
      !! H2O (one bond stretched)/6-31G CAS(4,4) SA-3, every pair
      type(error_type), allocatable, intent(out) :: error
      type(czt_molecule_t) :: mol
      type(casscf_result_t) :: result
      real(dp), parameter :: W(3) = [1.0_dp, 1.0_dp, 1.0_dp]/3.0_dp
      real(dp) :: d_csf(3, 3, 3), d_nocsf(3, 3, 3), h_csf(3, 3, 3)

      d_csf(:, :, 1) = H2O_12_D_CSF
      d_csf(:, :, 2) = H2O_13_D_CSF
      d_csf(:, :, 3) = H2O_23_D_CSF
      d_nocsf(:, :, 1) = H2O_12_D_NOCSF
      d_nocsf(:, :, 2) = H2O_13_D_NOCSF
      d_nocsf(:, :, 3) = H2O_23_D_NOCSF
      h_csf(:, :, 1) = H2O_12_H_CSF
      h_csf(:, :, 2) = H2O_13_H_CSF
      h_csf(:, :, 3) = H2O_23_H_CSF
      call converge([8, 1, 1], ["O ", "H ", "H "], H2O_XYZ, "6-31g", 10, 3, 4, 2, W, mol, &
                    result, error)
      if (allocated(error)) return
      call compare_pairs(mol, result, 3, 4, 2, W, reshape([1, 2, 1, 3, 2, 3], [2, 3]), &
                         d_csf, d_nocsf, h_csf, error)
      call mol%destroy()
   end subroutine test_cas44

end module test_mqc_sa_nac_pyscf_long

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_nac_pyscf_long, only: collect_mqc_sa_nac_pyscf_long_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_nac_pyscf_long", collect_mqc_sa_nac_pyscf_long_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
