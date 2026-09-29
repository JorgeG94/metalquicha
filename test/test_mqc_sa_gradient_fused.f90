!! The fused, all-roots-together SA-CASSCF gradient (phase 5)
module test_mqc_sa_gradient_fused
   !! Phase 5 of `SA_CASSCF_GRADIENT_PLAN.md`: `czt_sa_casscf_gradients`
   !! builds `build_sa_hessian` once and solves every requested root's
   !! Z-vector in one adaptive block PCG (`sa_block_zvector_solve`), instead
   !! of `czt_sa_casscf_gradient` called once per root. Gates:
   !!
   !! 1. `n_states = 1` is the same bit-identical short-circuit to
   !!    `czt_mcscf_gradient` the single-root entry point takes.
   !! 2. Per root, the fused path agrees with `czt_sa_casscf_gradient` (the
   !!    phase-4 single-root path) to <= 1e-10 max absolute component, on
   !!    LiH/STO-3G SA-2-CAS(2,2), LiH/6-31G CAS(4,4) SA-3 (all roots), and
   !!    C2H4 twisted 6-31G* SA-2 (`test_mqc_sa_gradient_pyscf_long.f90`'s
   !!    own geometry).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_physical_constants, only: BOHR_TO_ANGSTROM
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   use mqc_czt_mcscf_gradient, only: czt_mcscf_gradient
   use mqc_czt_sa_gradient, only: czt_sa_casscf_gradient, czt_sa_casscf_gradients
   implicit none
   private

   public :: collect_mqc_sa_gradient_fused_tests

   real(dp), parameter :: LIH(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 3.0139241961656_dp], [3, 2])
   integer, parameter :: LIH_Z(2) = [3, 1]
   character(len=2), parameter :: LIH_SYM(2) = ["Li", "H "]

   ! `test_mqc_sa_gradient_pyscf_long.f90`'s twisted C2H4, Angstrom.
   real(dp), parameter :: TWISTED_XYZ(3, 6) = reshape([ &
                                                      0.0_dp, 0.0_dp, 0.6695_dp, 0.0_dp, 0.0_dp, -0.6695_dp, &
                                                      0.0_dp, 0.9290_dp, 1.2320_dp, 0.0_dp, -0.9290_dp, 1.2320_dp, &
                                                   0.9290_dp, 0.1500_dp, -1.2320_dp, -0.9290_dp, 0.1500_dp, -1.2320_dp], &
                                                      [3, 6])
   integer, parameter :: C2H4_Z(6) = [6, 6, 1, 1, 1, 1]
   character(len=2), parameter :: C2H4_SYM(6) = ["C ", "C ", "H ", "H ", "H ", "H "]

contains

   subroutine collect_mqc_sa_gradient_fused_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("fused_n_states_1_is_bit_identical", test_short_circuit), &
                  new_unittest("fused_matches_single_root_lih_sto3g_sa2", test_lih_sa2), &
                  new_unittest("fused_matches_single_root_lih_631g_cas44_sa3", test_cas44_sa3), &
                  new_unittest("fused_matches_single_root_c2h4_twisted_sa2", test_c2h4_twisted) &
                  ]
   end subroutine collect_mqc_sa_gradient_fused_tests

   subroutine test_short_circuit(error)
      !! `czt_sa_casscf_gradients` at `n_states = 1` must equal
      !! `czt_mcscf_gradient` bit for bit, same as the single-root entry point
      !!
      !! At one thread: the densities are rebuilt from the CI vector, and a
      !! threaded RDM build's merge order is not fixed from call to call
      !! (`test_mqc_sa_gradient.f90`'s own bit-identity gate does the same).
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: result
      real(dp), allocatable :: g_base(:, :), g_fused(:, :, :)
      real(dp), allocatable :: ci3(:, :, :), en1(:)
      integer :: threads

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)

      call build_czt_molecule(LIH_Z, LIH_SYM, LIH, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "molecule should build")
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "RHF should converge")
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if
      call run_czt_casscf(mol, scf%orbitals, 1, 2, 1, 1, result, err, max_iterations=400, &
                          gradient_tol=1.0e-10_dp)
      call check(error,.not. err%has_error() .and. result%converged, &
                 "single-state CAS(2,2) should converge")
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if

      call czt_mcscf_gradient(mol, result%orbitals, 1, 2, result%dm1, result%dm2, g_base, err)
      call check(error,.not. err%has_error(), "czt_mcscf_gradient should not error")
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if

      allocate (ci3(size(result%ci_vector, 1), size(result%ci_vector, 2), 1))
      ci3(:, :, 1) = result%ci_vector
      en1 = [result%energy]
      call czt_sa_casscf_gradients(mol, result%orbitals, 1, 2, 1, 1, ci3, en1, [1.0_dp], [1], &
                                   g_fused, err)
      call check(error,.not. err%has_error(), "czt_sa_casscf_gradients should not error")
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if

      call omp_set_num_threads(threads)
      call check(error, all(g_fused(:, :, 1) == g_base), &
                 "the n_states=1 short-circuit should be bit-identical")
      call mol%destroy()
   end subroutine test_short_circuit

   subroutine test_lih_sa2(error)
      !! LiH/STO-3G SA-2-CAS(2,2): the fast, single-CI-direction-per-state case
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: WEIGHTS2(2) = [0.5_dp, 0.5_dp]
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: result
      real(dp), allocatable :: g_single(:, :), g_fused(:, :, :)
      integer :: root

      call build_czt_molecule(LIH_Z, LIH_SYM, LIH, "sto-3g", mol, err)
      call check(error,.not. err%has_error(), "molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "RHF should converge")
      if (allocated(error)) return
      call run_czt_casscf(mol, scf%orbitals, 1, 2, 1, 1, result, err, max_iterations=400, &
                          gradient_tol=1.0e-10_dp, n_states=2, weights=WEIGHTS2)
      call check(error,.not. err%has_error() .and. result%converged, "SA-2 should converge")
      if (allocated(error)) return

      call czt_sa_casscf_gradients(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                   result%energies, WEIGHTS2, [1, 2], g_fused, err)
      call check(error,.not. err%has_error(), "the fused gradient should build")
      if (allocated(error)) return

      do root = 1, 2
         call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                     result%energies, WEIGHTS2, root, g_single, err)
         call check(error,.not. err%has_error(), "the single-root gradient should build")
         if (allocated(error)) return
         call check(error, maxval(abs(g_fused(:, :, root) - g_single)) < 1.0e-10_dp, &
                    "fused and single-root gradients should agree to 1e-10")
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine test_lih_sa2

   subroutine test_cas44_sa3(error)
      !! LiH/6-31G CAS(4,4) SA-3: 18 non-redundant CI directions per state,
      !! all three roots requested from the fused entry point at once
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: WEIGHTS3(3) = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: result
      real(dp), allocatable :: g_single(:, :), g_fused(:, :, :)
      integer :: root

      call build_czt_molecule(LIH_Z, LIH_SYM, LIH, "6-31g", mol, err)
      call check(error,.not. err%has_error(), "molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "RHF should converge")
      if (allocated(error)) return
      call run_czt_casscf(mol, scf%orbitals, 0, 4, 2, 2, result, err, max_iterations=500, &
                          gradient_tol=1.0e-10_dp, n_states=3, weights=WEIGHTS3)
      call check(error,.not. err%has_error() .and. result%converged, "SA-3 should converge")
      if (allocated(error)) return

      call czt_sa_casscf_gradients(mol, result%orbitals, 0, 4, 2, 2, result%ci_vectors, &
                                   result%energies, WEIGHTS3, [1, 2, 3], g_fused, err)
      call check(error,.not. err%has_error(), "the fused gradient should build")
      if (allocated(error)) return

      do root = 1, 3
         call czt_sa_casscf_gradient(mol, result%orbitals, 0, 4, 2, 2, result%ci_vectors, &
                                     result%energies, WEIGHTS3, root, g_single, err)
         call check(error,.not. err%has_error(), "the single-root gradient should build")
         if (allocated(error)) return
         call check(error, maxval(abs(g_fused(:, :, root) - g_single)) < 1.0e-10_dp, &
                    "fused and single-root gradients should agree to 1e-10")
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine test_cas44_sa3

   subroutine test_c2h4_twisted(error)
      !! C2H4 twisted 6-31G* SA-2-CAS(2,2) -- a real polarised basis, and the
      !! symmetry-broken reference `developer_sa_casscf.rst` documents
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: WEIGHTS2(2) = [0.5_dp, 0.5_dp]
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      type(casscf_result_t) :: result
      real(dp), allocatable :: g_single(:, :), g_fused(:, :, :)
      integer :: root

      call build_czt_molecule(C2H4_Z, C2H4_SYM, TWISTED_XYZ/BOHR_TO_ANGSTROM, "6-31g_st_", &
                              mol, err)
      call check(error,.not. err%has_error(), "molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 16, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error,.not. err%has_error() .and. scf%converged, "RHF should converge")
      if (allocated(error)) return
      call run_czt_casscf(mol, scf%orbitals, 7, 2, 1, 1, result, err, max_iterations=400, &
                          gradient_tol=1.0e-10_dp, n_states=2, weights=WEIGHTS2)
      call check(error,.not. err%has_error() .and. result%converged, "SA-2 should converge")
      if (allocated(error)) return

      call czt_sa_casscf_gradients(mol, result%orbitals, 7, 2, 1, 1, result%ci_vectors, &
                                   result%energies, WEIGHTS2, [1, 2], g_fused, err)
      call check(error,.not. err%has_error(), "the fused gradient should build")
      if (allocated(error)) return

      do root = 1, 2
         call czt_sa_casscf_gradient(mol, result%orbitals, 7, 2, 1, 1, result%ci_vectors, &
                                     result%energies, WEIGHTS2, root, g_single, err)
         call check(error,.not. err%has_error(), "the single-root gradient should build")
         if (allocated(error)) return
         call check(error, maxval(abs(g_fused(:, :, root) - g_single)) < 1.0e-10_dp, &
                    "fused and single-root gradients should agree to 1e-10")
         if (allocated(error)) return
      end do
      call mol%destroy()
   end subroutine test_c2h4_twisted

end module test_mqc_sa_gradient_fused

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_gradient_fused, only: collect_mqc_sa_gradient_fused_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_gradient_fused", collect_mqc_sa_gradient_fused_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
