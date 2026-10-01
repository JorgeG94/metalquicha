!! Nonadiabatic couplings between SA-CASSCF states
module test_mqc_sa_nac
   !! LiH/STO-3G SA-2-CAS(2,2), the fast unit-test system shared with
   !! `test_mqc_sa_gradient`/`test_mqc_sa_gradient_pyscf_long`. Checks, in
   !! order:
   !!
   !! 1. `czt_sa_casscf_nac(1,2)` and `czt_sa_casscf_nac(2,1)` against PySCF
   !!    2.14's `pyscf.nac.sacasscf.NonAdiabaticCouplings`
   !!    (`tools/sa_casscf/pyscf_ref.py --nac`), both with and without the CSF
   !!    term, `d_IJ` and `h_IJ = (E_J-E_I) d_IJ`. The CI phase is arbitrary
   !!    (an overall sign per pair), so the comparison takes whichever of
   !!    `+PySCF`/`-PySCF` is closer and reports that sign.
   !! 2. `d_IJ = -d_JI` (antisymmetry) and `h_IJ` symmetric under `I<->J`.
   !! 3. Translational invariance of `h_IJ` **without** the CSF term (the CSF
   !!    term is not translationally invariant on its own; see the module
   !!    docstring on `mqc_czt_sa_nac`).
   !! 4. `czt_sa_casscf_nacs` (the shared-Hessian-state, several-pairs path)
   !!    reproduces `czt_sa_casscf_nac` (the single-pair path) for the same
   !!    pairs, to the Z-vector solver's tolerance.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   use mqc_czt_sa_nac, only: czt_sa_casscf_nac, czt_sa_casscf_nacs
   implicit none
   private

   public :: collect_mqc_sa_nac_tests

   real(dp), parameter :: WEIGHTS(2) = [0.5_dp, 0.5_dp]
   real(dp), parameter :: NAC_TOL = 5.0e-6_dp
      !! PySCF's own Z-vector solve is not tightened as far as the gradient
      !! gate's; loosened accordingly.
   real(dp), parameter :: LIH_XYZ(3, 2) = reshape([0.0_dp, 0.0_dp, 0.0_dp, &
                                                   0.0_dp, 0.0_dp, 1.5949_dp], [3, 2])

   !! PySCF reference, `tools/sa_casscf/pyscf_ref.py --system lih --basis
   !! sto-3g --nroots 2 --ncas 2 --nelecas 2 --nac 0,1` (0-based pair (0,1) =
   !! this module's (state_i=1, state_j=2)). Shape (3, natm): Li then H.
   !! `D_IJ` is PySCF's `kernel(mult_ediff=False)` verbatim -- a genuine
   !! physical NAC has no convention ambiguity beyond the CI phase. `H_IJ` is
   !! the **negative** of PySCF's raw `kernel(mult_ediff=True)`: PySCF scales
   !! by `e_bra - e_ket = E_I - E_J`, the opposite of this module's
   !! `h_IJ = (E_J - E_I) d_IJ` (`mqc_czt_sa_nac`'s module docstring), so its
   !! own number needs negating before it is this module's `h_IJ`.
   real(dp), parameter :: D12_WITH_CSF(3, 2) = reshape([ &
                                                       0.0_dp, 0.0_dp, 0.19332956788234_dp, &
                                                       0.0_dp, 0.0_dp, -0.07895003082021_dp], [3, 2])
   real(dp), parameter :: H12_WITH_CSF(3, 2) = reshape([ &
                                                       0.0_dp, 0.0_dp, 0.02496382071506_dp, &
                                                       0.0_dp, 0.0_dp, -0.01019448000858_dp], [3, 2])
   real(dp), parameter :: D12_NO_CSF(3, 2) = reshape([ &
                                                     0.0_dp, 0.0_dp, 0.09932700844621_dp, &
                                                     0.0_dp, 0.0_dp, -0.09932700844621_dp], [3, 2])
   real(dp), parameter :: H12_NO_CSF(3, 2) = reshape([ &
                                                     0.0_dp, 0.0_dp, 0.01282567202821_dp, &
                                                     0.0_dp, 0.0_dp, -0.01282567202821_dp], [3, 2])

contains

   subroutine collect_mqc_sa_nac_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lih_sto3g_sa2_nac_matches_pyscf", test_vs_pyscf), &
                  new_unittest("lih_sto3g_sa2_nac_identities", test_identities), &
                  new_unittest("lih_sto3g_sa2_nac_fused_matches_single", test_fused), &
                  new_unittest("degenerate_pair_is_refused", test_degenerate) &
                  ]
   end subroutine collect_mqc_sa_nac_tests

   subroutine converge_lih(mol, result, error)
      type(czt_molecule_t), intent(out) :: mol
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error

      type(rhf_result_t) :: scf

      call build_czt_molecule([3, 1], ["Li", "H "], LIH_XYZ/BOHR_TO_ANGSTROM, "sto-3g", &
                              mol, error)
      if (error%has_error()) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, error)
      if (error%has_error()) return
      call run_czt_casscf(mol, scf%orbitals, 1, 2, 1, 1, result, error, max_iterations=400, &
                          gradient_tol=1.0e-10_dp, n_states=2, weights=WEIGHTS)
   end subroutine converge_lih

   subroutine test_vs_pyscf(error)
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(casscf_result_t) :: result
      real(dp), allocatable :: coupling(:, :), interstate(:, :), csf(:, :)
      real(dp) :: ediff, sign_pick, worst

      call converge_lih(mol, result, err)
      call check(error,.not. err%has_error() .and. result%converged, &
                 "LiH SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call czt_sa_casscf_nac(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                             result%energies, WEIGHTS, 1, 2, coupling, interstate, csf, &
                             ediff, err, include_csf=.true.)
      call check(error,.not. err%has_error(), "the (1,2) NAC with CSF should build")
      if (allocated(error)) return

      ! The CI phase is arbitrary: take whichever overall sign agrees.
      sign_pick = merge(1.0_dp, -1.0_dp, &
                        sum(coupling*D12_WITH_CSF) >= 0.0_dp)
      worst = maxval(abs(sign_pick*coupling - D12_WITH_CSF))
      write (*, "(a,es10.2)") "    d_12 (with CSF) max |mqc - PySCF| = ", worst
      call check(error, worst < NAC_TOL, "d_12 with CSF should match PySCF")
      if (allocated(error)) return
      worst = maxval(abs(sign_pick*interstate - H12_WITH_CSF))
      write (*, "(a,es10.2)") "    h_12 (with CSF) max |mqc - PySCF| = ", worst
      call check(error, worst < NAC_TOL, "h_12 with CSF should match PySCF")
      if (allocated(error)) return

      call czt_sa_casscf_nac(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                             result%energies, WEIGHTS, 1, 2, coupling, interstate, csf, &
                             ediff, err, include_csf=.false.)
      call check(error,.not. err%has_error(), "the (1,2) NAC without CSF should build")
      if (allocated(error)) return
      worst = maxval(abs(sign_pick*coupling - D12_NO_CSF))
      write (*, "(a,es10.2)") "    d_12 (no CSF) max |mqc - PySCF| = ", worst
      call check(error, worst < NAC_TOL, "d_12 without CSF should match PySCF")
      if (allocated(error)) return
      worst = maxval(abs(sign_pick*interstate - H12_NO_CSF))
      write (*, "(a,es10.2)") "    h_12 (no CSF) max |mqc - PySCF| = ", worst
      call check(error, worst < NAC_TOL, "h_12 without CSF should match PySCF")
      if (allocated(error)) return

      call mol%destroy()
   end subroutine test_vs_pyscf

   subroutine test_identities(error)
      !! `d_IJ = -d_JI`, `h_IJ` symmetric, and translational invariance of
      !! `h_IJ` without the CSF term
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(casscf_result_t) :: result
      real(dp), allocatable :: coupling12(:, :), interstate12(:, :), csf12(:, :)
      real(dp), allocatable :: coupling21(:, :), interstate21(:, :), csf21(:, :)
      real(dp) :: e12, e21, worst

      call converge_lih(mol, result, err)
      call check(error,.not. err%has_error() .and. result%converged, &
                 "LiH SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call czt_sa_casscf_nac(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                             result%energies, WEIGHTS, 1, 2, coupling12, interstate12, csf12, &
                             e12, err)
      call check(error,.not. err%has_error(), "the (1,2) NAC should build")
      if (allocated(error)) return
      call czt_sa_casscf_nac(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                             result%energies, WEIGHTS, 2, 1, coupling21, interstate21, csf21, &
                             e21, err)
      call check(error,.not. err%has_error(), "the (2,1) NAC should build")
      if (allocated(error)) return

      call check(error, abs(e12 + e21) < 1.0e-12_dp, "E_J - E_I should flip sign")
      if (allocated(error)) return

      worst = maxval(abs(coupling12 + coupling21))
      write (*, "(a,es10.2)") "    |d_12 + d_21| (should be ~0) = ", worst
      call check(error, worst < NAC_TOL, "d_IJ should be antisymmetric")
      if (allocated(error)) return

      worst = maxval(abs(interstate12 - interstate21))
      write (*, "(a,es10.2)") "    |h_12 - h_21| (should be ~0) = ", worst
      call check(error, worst < NAC_TOL, "h_IJ should be symmetric")
      if (allocated(error)) return

      worst = maxval(abs(sum(interstate12 - csf12, dim=2)))
      write (*, "(a,es10.2)") "    |sum_atoms (h_12 - csf_12)| (should be ~0) = ", worst
      call check(error, worst < NAC_TOL, &
                 "h_IJ without the CSF term should be translationally invariant")
      if (allocated(error)) return

      call mol%destroy()
   end subroutine test_identities

   subroutine test_degenerate(error)
      !! A pair whose energies coincide is refused by name, not divided by zero
      !!
      !! The energies are made degenerate by hand on a real converged state:
      !! nothing before the division reads them except `build_sa_hessian`'s
      !! active energies, which this refusal precedes.
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(casscf_result_t) :: result
      real(dp), allocatable :: coupling(:, :), interstate(:, :), csf(:, :)
      real(dp), allocatable :: energies(:)
      real(dp) :: ediff

      call converge_lih(mol, result, err)
      call check(error,.not. err%has_error() .and. result%converged, &
                 "LiH SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      energies = result%energies
      energies(2) = energies(1)
      call czt_sa_casscf_nac(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                             energies, WEIGHTS, 1, 2, coupling, interstate, csf, ediff, err)
      call check(error, err%has_error(), "a degenerate pair should be refused")
      if (allocated(error)) return
      call check(error, err%get_code() == ERROR_VALIDATION, &
                 "and refused as a validation error: "//err%get_message())
      if (allocated(error)) return
      if (allocated(coupling)) then
         call check(error, all(coupling == coupling), "no NaN should be left in the coupling")
      end if
   end subroutine test_degenerate

   subroutine test_fused(error)
      !! `czt_sa_casscf_nacs`, sharing one SA Hessian state across both pairs,
      !! reproduces `czt_sa_casscf_nac` called once per pair
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(casscf_result_t) :: result
      real(dp), allocatable :: coupling(:, :), interstate(:, :), csf(:, :)
      real(dp), allocatable :: couplings(:, :, :), interstates(:, :, :), csfs(:, :, :)
      real(dp), allocatable :: ediffs(:)
      integer, allocatable :: pairs(:, :)
      real(dp) :: ediff, worst

      call converge_lih(mol, result, err)
      call check(error,.not. err%has_error() .and. result%converged, &
                 "LiH SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      allocate (pairs(2, 2))
      pairs(:, 1) = [1, 2]
      pairs(:, 2) = [2, 1]
      call czt_sa_casscf_nacs(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                              result%energies, WEIGHTS, pairs, couplings, interstates, csfs, &
                              ediffs, err)
      call check(error,.not. err%has_error(), "the fused pair list should build")
      if (allocated(error)) return

      call czt_sa_casscf_nac(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                             result%energies, WEIGHTS, 1, 2, coupling, interstate, csf, &
                             ediff, err)
      call check(error,.not. err%has_error(), "the single-pair (1,2) NAC should build")
      if (allocated(error)) return

      worst = maxval(abs(couplings(:, :, 1) - coupling))
      write (*, "(a,es10.2)") "    fused vs single, d_12 max diff = ", worst
      call check(error, worst < 1.0e-9_dp, "the fused path should match the single-pair one")
      if (allocated(error)) return
      call check(error, abs(ediffs(1) - ediff) < 1.0e-12_dp, "energy gaps should match")
      if (allocated(error)) return

      call mol%destroy()
   end subroutine test_fused

end module test_mqc_sa_nac

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_nac, only: collect_mqc_sa_nac_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_nac", collect_mqc_sa_nac_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
