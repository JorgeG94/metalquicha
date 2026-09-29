!! State-averaged CASSCF: energy, singlet selection, and the n_states=1 seam
module test_mqc_sa_casscf
   !! Phase 1 of the SA-CASSCF project (`SA_CASSCF_GRADIENT_PLAN.md` at the
   !! repository root; equations in `mqc_docs/source/developer_sa_casscf.rst`).
   !!
   !! **`n_states = 1` has to be bit-identical to `run_czt_casscf` without the
   !! optional arguments at all.** State averaging shares that routine's macro
   !! loop and branches only where the CI is solved, the densities are built
   !! and the objective is read, so a single state takes exactly the old
   !! operations. `sa_n_states_one_matches_plain_path` checks this directly.
   !!
   !! References (LiH, STO-3G, CAS(2,2), Li at the origin, H at
   !! `z = 1.5949` Angstrom): PySCF 2.14, `mcscf.CASSCF(mf, 2,
   !! 2).state_average_([0.5, 0.5])`, fed this repository's own basis JSON --
   !! see `tools/sa_casscf/pyscf_ref.py`. `SA_CASSCF_GRADIENT_PLAN.md`'s
   !! progress log carries the same two numbers.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_casci, only: run_czt_casci, casci_result_t
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t, generalized_fock, &
                            orbital_gradient, mcscf_fock_t
   use mqc_rdm, only: spin_squared
   implicit none
   private

   public :: collect_mqc_sa_casscf_tests

   real(dp), parameter :: LIH(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 3.0139241961656_dp], [3, 2])
      !! `1.5949` Angstrom in Bohr (`/0.529177210903`), matching
      !! `tools/sa_casscf/lih.xyz`.
   integer, parameter :: LIH_Z(2) = [3, 1]
   character(len=2), parameter :: LIH_SYM(2) = ["Li", "H "]

   real(dp), parameter :: E_ROOT1 = -7.854944206511_dp
   real(dp), parameter :: E_ROOT2 = -7.725818481002_dp

contains

   subroutine collect_mqc_sa_casscf_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("lih_sa2_energy_and_spin", test_lih_sa2), &
                  new_unittest("n_states_one_matches_plain_path", &
                               test_n_states_one_matches_plain), &
                  new_unittest("singlet_selection_excludes_the_triplet", &
                               test_singlet_selection), &
                  new_unittest("sa_orbital_gradient_vanishes_root_gradient_does_not", &
                               test_sa_gradient_not_root_gradient) &
                  ]
   end subroutine collect_mqc_sa_casscf_tests

   subroutine lih_reference(mol, orbitals, err)
      !! LiH/STO-3G, converged RHF orbitals
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      type(error_t), intent(inout) :: err

      type(rhf_result_t) :: scf

      call build_czt_molecule(LIH_Z, LIH_SYM, LIH, "sto-3g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error()) return
      if (.not. scf%converged) then
         call err%set(ERROR_VALIDATION, "the LiH RHF reference did not converge")
         return
      end if
      orbitals = scf%orbitals
   end subroutine lih_reference

   subroutine test_lih_sa2(error)
      !! SA-2-CAS(2,2) against PySCF, and every root a singlet
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), parameter :: WEIGHTS(2) = [0.5_dp, 0.5_dp]

      call lih_reference(mol, orbitals, err)
      call check(error,.not. err%has_error(), "the reference SCF should converge")
      if (allocated(error)) return

      call run_czt_casscf(mol, orbitals, 1, 2, 1, 1, result, err, &
                          max_iterations=300, gradient_tol=1.0e-10_dp, &
                          n_states=2, weights=WEIGHTS)
      call mol%destroy()
      call check(error,.not. err%has_error(), "SA-CASSCF should not error")
      if (allocated(error)) return
      call check(error, result%converged, "the orbital optimisation should converge")
      if (allocated(error)) return

      call check(error, allocated(result%energies), "every root's energy is carried")
      if (allocated(error)) return
      call check(error, size(result%energies), 2, "two states were asked for")
      if (allocated(error)) return

      call check(error, result%energies(1), E_ROOT1, "root 1 vs PySCF", thr=1.0e-8_dp)
      if (allocated(error)) return
      call check(error, result%energies(2), E_ROOT2, "root 2 vs PySCF", thr=1.0e-8_dp)
      if (allocated(error)) return
      call check(error, result%energy, 0.5_dp*(E_ROOT1 + E_ROOT2), &
                 "E_SA is the equal-weight average", thr=1.0e-8_dp)
      if (allocated(error)) return

      call check(error, allocated(result%spins), "every root's <S^2> is carried")
      if (allocated(error)) return
      call check(error, result%spins(1), 0.0_dp, "root 1 is a singlet", thr=1.0e-8_dp)
      if (allocated(error)) return
      call check(error, result%spins(2), 0.0_dp, "root 2 is a singlet", thr=1.0e-8_dp)
   end subroutine test_lih_sa2

   subroutine test_n_states_one_matches_plain(error)
      !! `n_states = 1, weights = [1]` reproduces the no-optional-argument
      !! path bit for bit
      !!
      !! `run_czt_casscf` branches on `n_states > 1` only where the CI is
      !! solved, the densities are built and the objective is read. This test
      !! checks that `n_states = 1` takes the single-state operations.
      !!
      !! The threshold is `1e-9`, not exactly zero: two separate calls to the
      !! *same* code path already scatter by ~5e-10 under OpenMP, because a
      !! threaded reduction's merge order is not fixed run to run (see this
      !! repository's own note on that). That noise has nothing to do with
      !! `n_states`, so asking for literal bit-identity here would be testing
      !! the thread count rather than the dispatch. `1e-9` is ten orders of
      !! magnitude tighter than the `1e-8` a converged energy is judged by
      !! elsewhere, which is as close to "no bit changes" as a threaded build
      !! can promise.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: plain, with_n_states
      real(dp), parameter :: ONE_WEIGHT(1) = [1.0_dp]

      call lih_reference(mol, orbitals, err)
      call check(error,.not. err%has_error(), "the reference SCF should converge")
      if (allocated(error)) return

      call run_czt_casscf(mol, orbitals, 1, 2, 1, 1, plain, err, &
                          max_iterations=300, gradient_tol=1.0e-10_dp)
      call check(error,.not. err%has_error(), "the plain path should not error")
      if (allocated(error)) return

      call run_czt_casscf(mol, orbitals, 1, 2, 1, 1, with_n_states, err, &
                          max_iterations=300, gradient_tol=1.0e-10_dp, &
                          n_states=1, weights=ONE_WEIGHT)
      call mol%destroy()
      call check(error,.not. err%has_error(), "n_states=1 should not error")
      if (allocated(error)) return

      call check(error, with_n_states%energy, plain%energy, &
                 "n_states=1 reproduces the plain path", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, with_n_states%gradient_norm, plain%gradient_norm, &
                 "n_states=1 reproduces the plain path's gradient norm", thr=1.0e-9_dp)
      if (allocated(error)) return
      call check(error, with_n_states%iterations, plain%iterations, &
                 "n_states=1 takes the same macro-iterations")
      if (allocated(error)) return
      call check(error, maxval(abs(with_n_states%orbitals - plain%orbitals)), 0.0_dp, &
                 "n_states=1 reproduces the plain path's orbitals", thr=1.0e-9_dp)
   end subroutine test_n_states_one_matches_plain

   subroutine test_singlet_selection(error)
      !! Every determinant in a CAS(2,2) with one active electron per spin
      !! decomposes into three singlets and one triplet (Sec. "Spin",
      !! `mqc_docs/source/developer_sa_casscf.rst`). Unsymmetrised, all four
      !! roots are reachable and one of them is the triplet; symmetrised, the
      !! Davidson never sees it.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casci_result_t) :: unsymmetrized, symmetrized
      real(dp) :: s2
      integer :: i
      logical :: found_triplet

      call lih_reference(mol, orbitals, err)
      call check(error,.not. err%has_error(), "the reference SCF should converge")
      if (allocated(error)) return

      ! All four determinants of this active space, unrestricted: three
      ! singlets and the Ms=0 triplet component, per the note above.
      call run_czt_casci(mol, orbitals, 1, 2, 1, 1, unsymmetrized, err, n_roots=4, &
                         tolerance=1.0e-11_dp)
      call check(error,.not. err%has_error(), "the unsymmetrised CASCI should not error")
      if (allocated(error)) return

      found_triplet = .false.
      do i = 1, 4
         s2 = spin_squared(2, 1, 1, unsymmetrized%vectors(:, :, i), err)
         call check(error,.not. err%has_error(), "spin_squared should not error")
         if (allocated(error)) return
         if (abs(s2 - 2.0_dp) < 1.0e-6_dp) found_triplet = .true.
      end do
      call check(error, found_triplet, &
                 "one of the four unsymmetrised roots should be the triplet")
      if (allocated(error)) return

      ! Symmetrised, two roots: the two lowest singlets, both with <S^2> = 0.
      call run_czt_casci(mol, orbitals, 1, 2, 1, 1, symmetrized, err, n_roots=2, &
                         tolerance=1.0e-11_dp, symmetrize_singlet=.true.)
      call mol%destroy()
      call check(error,.not. err%has_error(), "the symmetrised CASCI should not error")
      if (allocated(error)) return

      do i = 1, 2
         s2 = spin_squared(2, 1, 1, symmetrized%vectors(:, :, i), err)
         call check(error,.not. err%has_error(), "spin_squared should not error")
         if (allocated(error)) return
         call check(error, s2, 0.0_dp, "a symmetrised root should be a singlet", &
                    thr=1.0e-8_dp)
         if (allocated(error)) return
      end do
   end subroutine test_singlet_selection

   subroutine test_sa_gradient_not_root_gradient(error)
      !! At SA convergence, `orbital_gradient` on the SA densities is ~0; on
      !! root 1's own densities it is not -- the sanity check that the
      !! optimiser drove `E_SA` to a stationary point and not `E_1`.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: sa_result
      type(casci_result_t) :: root_ci
      type(mcscf_fock_t) :: fock_sa, fock_root1
      real(dp), allocatable :: gradient_sa(:, :), gradient_root1(:, :)
      real(dp), parameter :: WEIGHTS(2) = [0.5_dp, 0.5_dp]

      call lih_reference(mol, orbitals, err)
      call check(error,.not. err%has_error(), "the reference SCF should converge")
      if (allocated(error)) return

      call run_czt_casscf(mol, orbitals, 1, 2, 1, 1, sa_result, err, &
                          max_iterations=300, gradient_tol=1.0e-10_dp, &
                          n_states=2, weights=WEIGHTS)
      call check(error,.not. err%has_error(), "SA-CASSCF should not error")
      if (allocated(error)) return
      call check(error, sa_result%converged, "the SA optimisation should converge")
      if (allocated(error)) return

      ! The SA densities the optimiser actually used -- `orbital_gradient` on
      ! these is what `sa_result%gradient_norm` already reports as converged.
      call generalized_fock(mol, sa_result%orbitals, 1, 2, sa_result%dm1, &
                            sa_result%dm2, fock_sa, err)
      call check(error,.not. err%has_error(), "generalized_fock (SA) should not error")
      if (allocated(error)) return
      call orbital_gradient(fock_sa, 1, 2, gradient_sa)
      call check(error, maxval(abs(gradient_sa)) < 1.0e-8_dp, &
                 "the SA orbital gradient is ~0 at the SA-converged orbitals")
      if (allocated(error)) return

      ! Root 1's own density, at the SAME converged orbitals: a CASCI, not a
      ! CASSCF, so its own orbital gradient need not vanish there.
      call run_czt_casci(mol, sa_result%orbitals, 1, 2, 1, 1, root_ci, err, &
                         n_roots=2, tolerance=1.0e-11_dp, symmetrize_singlet=.true.)
      call check(error,.not. err%has_error(), "the root-1 CASCI should not error")
      if (allocated(error)) return

      call generalized_fock(mol, sa_result%orbitals, 1, 2, root_ci%dm1, root_ci%dm2, &
                            fock_root1, err)
      call mol%destroy()
      call check(error,.not. err%has_error(), "generalized_fock (root 1) should not error")
      if (allocated(error)) return
      call orbital_gradient(fock_root1, 1, 2, gradient_root1)
      call check(error, maxval(abs(gradient_root1)) > 1.0e-4_dp, &
                 "root 1's own orbital gradient is not zero at the SA optimum")
   end subroutine test_sa_gradient_not_root_gradient

end module test_mqc_sa_casscf

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_casscf, only: collect_mqc_sa_casscf_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_casscf", collect_mqc_sa_casscf_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
