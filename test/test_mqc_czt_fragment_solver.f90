!! The fragment-solver interface, against the backend calls it replaces
module test_mqc_czt_fragment_solver
   !! `solve_fragment_method` collapses FMO's, EE-MBE's and EFMO's own
   !! `run_czt_rhf`/`run_czt_mp2` calls into one interface; see
   !! `mqc_docs/source/developer_fragment_solver.rst`. What is checked here is
   !! that it reproduces those direct calls, on water in STO-3G, fast:
   !!
   !!   1. Plain Hartree-Fock through the solver matches a direct `run_czt_rhf`
   !!      call, energy and density both.
   !!   2. With an embedding operator, `internal` is `energy - Tr(D u)` *exactly*
   !!      -- one run's own arithmetic -- and `energy` still matches a direct
   !!      `run_czt_rhf(..., h_extra=)` call.
   !!   3. MP2 from an already-converged `reference` runs no SCF -- `outcome%scf`
   !!      is that reference, exactly, with no arithmetic done to it -- and
   !!      `correlation` matches a direct `run_czt_mp2` call on its orbitals,
   !!      same-spin plus opposite-spin.
   !!   4. A functional is refused, by an internal-consistency message, rather
   !!      than answered as Hartree-Fock.
   !!
   !! **"Matches" rather than "to the bit" wherever two SCFs or two MP2s are
   !! run separately** -- the solver's own and a direct call made again here to
   !! check it -- and compared with `SCF_TOL` rather than `thr=0.0`. Each is
   !! deterministic on its own, but this repository's own bit-identity
   !! guarantee (`mqc_docs/source/developer_fragment_solver.rst`) is stated at
   !! one thread; this test suite runs at up to four, where a threaded BLAS
   !! reduction need not reassociate identically between two independent runs.
   !! A same-run scalar expression on values one call already returned has no
   !! such freedom and is asserted exactly.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use mqc_error, only: error_t
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mp2, only: mp2_result_t, run_czt_mp2
   use mqc_czt_esp, only: esp_matrices
   use mqc_cuest_iface, only: cuest_scf_settings_t
   use mqc_czt_fragment_solver, only: fragment_request_t, fragment_outcome_t, &
                                      solve_fragment_method
   use mqc_physical_constants, only: ANGSTROM_TO_BOHR
   implicit none
   private

   public :: collect_mqc_czt_fragment_solver_tests

   real(dp), parameter :: ANG = ANGSTROM_TO_BOHR
   character(len=*), parameter :: BASIS = "sto-3g"
   integer, parameter :: NELEC = 10   !! Water: 8 + 1 + 1
   real(dp), parameter :: SCF_TOL = 1.0e-11_dp
      !! For comparing two SEPARATELY run SCFs -- the solver's and a direct
      !! `run_czt_rhf` call. Each is deterministic on its own, but a threaded
      !! BLAS reduction need not reassociate the same way between two runs, so
      !! "to the bit" is asserted only within one run's own arithmetic (a
      !! scalar expression on values that call already returned); across two
      !! independent SCFs this is what stands in for it.

contains

   subroutine collect_mqc_czt_fragment_solver_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("solver_plain_hf_matches_run_czt_rhf_to_the_bit", &
                               test_plain_hf), &
                  new_unittest("solver_embedded_hf_internal_is_energy_minus_tr_d_u", &
                               test_embedded_hf), &
                  new_unittest("solver_mp2_from_a_reference_runs_no_scf", &
                               test_mp2_reference), &
                  new_unittest("solver_refuses_a_functional", test_refuses_functional) &
                  ]
   end subroutine collect_mqc_czt_fragment_solver_tests

   subroutine water_geometry(z, symbols, xyz)
      !! One water, Bohr, in the `yz` plane
      integer, intent(out) :: z(3)
      character(len=2), intent(out) :: symbols(3)
      real(dp), intent(out) :: xyz(3, 3)

      z = [8, 1, 1]
      symbols = ["O ", "H ", "H "]
      xyz = reshape([0.00000000_dp, 0.00000000_dp, 0.11726921_dp, &
                     0.00000000_dp, 0.75698224_dp, -0.46907684_dp, &
                     0.00000000_dp, -0.75698224_dp, -0.46907684_dp], [3, 3])*ANG
   end subroutine water_geometry

   subroutine test_plain_hf(error)
      !! Plain Hartree-Fock through the solver equals a direct `run_czt_rhf`
      type(error_type), allocatable, intent(out) :: error

      integer :: z(3)
      character(len=2) :: symbols(3)
      real(dp) :: xyz(3, 3)
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(cuest_scf_settings_t) :: method
      type(fragment_request_t) :: request
      type(fragment_outcome_t) :: outcome
      type(rhf_result_t) :: direct

      call water_geometry(z, symbols, xyz)
      call build_czt_molecule(z, symbols, xyz, BASIS, mol, err)
      call check(error,.not. err%has_error(), "building the molecule failed: "// &
                 err%get_message())
      if (allocated(error)) return

      request%max_iter = 100
      request%energy_tol = 1.0e-10_dp
      request%density_tol = 1.0e-8_dp

      call solve_fragment_method(method, mol, NELEC, z, request, outcome, err)
      call check(error,.not. err%has_error(), "solve_fragment_method failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, outcome%converged, "the solver's SCF did not converge")
      if (allocated(error)) return

      ! Exactly the call `solve_fragment_method` makes: every optional
      ! unallocated, so it arrives absent at `run_czt_rhf` on both sides.
      call run_czt_rhf(mol, NELEC, request%max_iter, request%energy_tol, &
                       request%density_tol, request%verbose, direct, err, &
                       scf=request%drive, guess=request%guess, &
                       guess_density=request%guess_density, h_extra=request%h_extra, &
                       projector=request%projector, grad_tol=request%grad_tol)
      call check(error,.not. err%has_error(), "the direct RHF failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call check(error, outcome%energy, direct%energy, thr=SCF_TOL, &
                 message="the solver's energy is not run_czt_rhf's")
      if (allocated(error)) return
      ! Within this one run's own arithmetic, so exact: `reference` and
      ! `internal` are `energy` verbatim with no `h_extra` to subtract.
      call check(error, outcome%reference, outcome%energy, thr=0.0_dp, &
                 message="reference should equal energy on a plain Hartree-Fock request")
      if (allocated(error)) return
      call check(error, outcome%internal, outcome%energy, thr=0.0_dp, &
                 message="internal should equal energy with no h_extra")
      if (allocated(error)) return
      call check(error, maxval(abs(outcome%density - direct%density)) < SCF_TOL, &
                 "the solver's density is not run_czt_rhf's")
      if (allocated(error)) return
      call mol%destroy()
   end subroutine test_plain_hf

   subroutine test_embedded_hf(error)
      !! With an embedding operator, `internal` is `energy - Tr(D u)`, and
      !! `energy` still matches a direct `run_czt_rhf(..., h_extra=)`
      type(error_type), allocatable, intent(out) :: error

      integer :: z(3)
      character(len=2) :: symbols(3)
      real(dp) :: xyz(3, 3), points(3, 1)
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(cuest_scf_settings_t) :: method
      type(fragment_request_t) :: request
      type(fragment_outcome_t) :: outcome
      type(rhf_result_t) :: direct
      real(dp), allocatable :: matrices(:, :, :)

      call water_geometry(z, symbols, xyz)
      call build_czt_molecule(z, symbols, xyz, BASIS, mol, err)
      call check(error,.not. err%has_error(), "building the molecule failed: "// &
                 err%get_message())
      if (allocated(error)) return

      ! A point charge off to the side: any symmetric one-electron operator
      ! over the molecule's own AOs does, and this is what the callers build
      ! theirs from.
      points(:, 1) = [3.0_dp, 0.0_dp, 0.0_dp]*ANG
      call esp_matrices(mol, points, matrices, err)
      call check(error,.not. err%has_error(), "building h_extra failed: "// &
                 err%get_message())
      if (allocated(error)) return

      request%max_iter = 100
      request%energy_tol = 1.0e-10_dp
      request%density_tol = 1.0e-8_dp
      request%h_extra = matrices(:, :, 1)

      call solve_fragment_method(method, mol, NELEC, z, request, outcome, err)
      call check(error,.not. err%has_error(), "solve_fragment_method failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, outcome%converged, "the solver's SCF did not converge")
      if (allocated(error)) return

      call check(error, outcome%internal, &
                 outcome%energy - sum(outcome%density*request%h_extra), thr=0.0_dp, &
                 message="internal is not energy - Tr(D u), to the bit")
      if (allocated(error)) return

      call run_czt_rhf(mol, NELEC, request%max_iter, request%energy_tol, &
                       request%density_tol, request%verbose, direct, err, &
                       scf=request%drive, guess=request%guess, &
                       guess_density=request%guess_density, h_extra=request%h_extra, &
                       projector=request%projector, grad_tol=request%grad_tol)
      call check(error,.not. err%has_error(), "the direct RHF failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, outcome%energy, direct%energy, thr=SCF_TOL, &
                 message="the embedded solver energy is not run_czt_rhf's")
      if (allocated(error)) return
      call mol%destroy()
   end subroutine test_embedded_hf

   subroutine test_mp2_reference(error)
      !! MP2 from an already-converged `reference` runs no SCF
      type(error_type), allocatable, intent(out) :: error

      integer :: z(3)
      character(len=2) :: symbols(3)
      real(dp) :: xyz(3, 3)
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(cuest_scf_settings_t) :: method
      type(fragment_request_t) :: request
      type(fragment_outcome_t) :: outcome
      type(rhf_result_t) :: reference
      type(mp2_result_t) :: mp2_direct

      call water_geometry(z, symbols, xyz)
      call build_czt_molecule(z, symbols, xyz, BASIS, mol, err)
      call check(error,.not. err%has_error(), "building the molecule failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call run_czt_rhf(mol, NELEC, 100, 1.0e-10_dp, 1.0e-8_dp, .false., reference, err)
      call check(error,.not. err%has_error(), "the reference RHF failed: "// &
                 err%get_message())
      if (allocated(error)) return
      call check(error, reference%converged, "the reference SCF did not converge")
      if (allocated(error)) return

      ! No core to freeze, so the direct `run_czt_mp2` below needs no
      ! `n_frozen` either: both correlate every occupied orbital.
      method%run_mp2 = .true.
      method%freeze_core = .false.
      call solve_fragment_method(method, mol, NELEC, z, request, outcome, err, &
                                 reference=reference)
      call check(error,.not. err%has_error(), "solve_fragment_method failed: "// &
                 err%get_message())
      if (allocated(error)) return

      ! `outcome%scf` is `reference`, unchanged: no SCF ran here.
      call check(error, outcome%scf%energy, reference%energy, thr=0.0_dp, &
                 message="a present reference should skip the SCF")
      if (allocated(error)) return
      call check(error, outcome%reference, reference%energy, thr=0.0_dp, &
                 message="outcome%reference is not the reference's energy")
      if (allocated(error)) return

      call run_czt_mp2(mol, reference%orbitals, reference%orbital_energies, NELEC/2, &
                       reference%energy, mp2_direct, err)
      call check(error,.not. err%has_error(), "the direct MP2 failed: "// &
                 err%get_message())
      if (allocated(error)) return

      call check(error, outcome%correlation, mp2_direct%same_spin + mp2_direct%opposite_spin, &
                 thr=SCF_TOL, &
                 message="the solver's correlation is not run_czt_mp2's")
      if (allocated(error)) return
      call check(error, outcome%energy, reference%energy + outcome%correlation, thr=0.0_dp, &
                 message="energy is not reference + correlation")
      if (allocated(error)) return
      call mol%destroy()
   end subroutine test_mp2_reference

   subroutine test_refuses_functional(error)
      !! A functional request is an internal-consistency error, not DFT
      type(error_type), allocatable, intent(out) :: error

      integer :: z(3)
      type(czt_molecule_t) :: mol
      type(error_t) :: err
      type(cuest_scf_settings_t) :: method
      type(fragment_request_t) :: request
      type(fragment_outcome_t) :: outcome
      character(len=:), allocatable :: message

      z = [8, 1, 1]
      method%functional = "pbe"

      ! `mol` is never read on this path -- the refusal comes before it would
      ! be -- so it is passed unbuilt.
      call solve_fragment_method(method, mol, NELEC, z, request, outcome, err)
      call check(error, err%has_error(), &
                 "a functional request should be refused rather than run as Hartree-Fock")
      if (allocated(error)) return

      message = err%get_message()
      call check(error, index(message, "fragment_refusal") > 0, &
                 "the refusal should point at fragment_refusal: "//message)
   end subroutine test_refuses_functional

end module test_mqc_czt_fragment_solver

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_czt_fragment_solver, only: collect_mqc_czt_fragment_solver_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_czt_fragment_solver", &
                               collect_mqc_czt_fragment_solver_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
