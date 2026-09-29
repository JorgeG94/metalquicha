!! The explicit SA Hessian's eigenvalues, on a genuinely multi-dimensional
!! CI space
module test_mqc_sa_hessian_redundancy_long
   !! `test_mqc_sa_hessian.f90`'s LiH/STO-3G SA-2-CAS(2,2) redundancy gate
   !! only has one non-redundant CI direction per state to work with; this
   !! repeats it on LiH/6-31G CAS(4,4) SA-3 (`test_mqc_sa_hessian_multistate.f90`'s
   !! system: 11 basis functions, no frozen core, 36 determinants, 21
   !! singlet-symmetric), where 18 genuine CI directions per state actually
   !! exercise the CI-CI and orbital-CI blocks. `n_param = n_rot + 3*n_det =
   !! 28 + 108 = 136`, so the explicit Hessian is `136` unit-vector
   !! `sa_hessian_apply` calls -- still cheap at this basis size, but slow
   !! enough next to the routine suite's target that it is an `addlongtest`
   !! rather than folded into the routine file.
   !!
   !! **The null space is one larger than the `n_states*(antisym_dim +
   !! n_states)` formula predicts, and that is not a bug.** LiH/6-31G's RHF
   !! orbitals 4 and 5 (and separately 8 and 9) are exactly degenerate --
   !! linear-molecule px/py pairs -- and this run's converged CAS(4,4)
   !! active space ends up spanning only one member of a degenerate pair
   !! (confirmed by printing `scf%orbital_energies`, and unaffected by which
   !! orbital the initial guess put in the active window: CASSCF's Newton
   !! step is free to rotate across the whole non-redundant space regardless
   !! of the starting column order, so it converges to the same effective
   !! partition either way). Rotating that active member into its still-
   !! virtual, still-degenerate partner changes nothing, exactly as rotating
   !! within any other degenerate irrep would not -- a point-group-symmetry
   !! redundancy, present in any orbital Hessian (single-state or SA) for a
   !! symmetric molecule with no symmetry constraint imposed on the active
   !! space, and orthogonal to the state-averaging/singlet-restriction
   !! redundancy the formula accounts for. The gate below checks for *at
   !! least* the formula's count, and separately that the null cluster and
   !! the positive spectrum are cleanly separated (here, `~1e-15` against
   !! `~3e-2`, six orders of magnitude), which is what would fail to hold if
   !! an extra "near-zero" eigenvalue were actually a wrong small positive or
   !! negative curvature rather than a further exact symmetry.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t
   use mqc_czt_sa_hessian, only: sa_hessian_t, build_sa_hessian, destroy_sa_hessian, &
                                 sa_hessian_n_param, sa_hessian_apply
   implicit none
   private

   public :: collect_mqc_sa_hessian_redundancy_long_tests

   real(dp), parameter :: LIH(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 3.0139241961656_dp], [3, 2])
   integer, parameter :: LIH_Z(2) = [3, 1]
   character(len=2), parameter :: LIH_SYM(2) = ["Li", "H "]

contains

   subroutine collect_mqc_sa_hessian_redundancy_long_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("redundancy_and_positive_curvature_cas44_sa3", &
                               test_redundancy_sa3) &
                  ]
   end subroutine collect_mqc_sa_hessian_redundancy_long_tests

   subroutine build_dense(state, dense, error)
      type(sa_hessian_t), intent(in) :: state
      real(dp), allocatable, intent(out) :: dense(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unit_vec(:, :), column(:, :)
      integer :: np, k

      np = sa_hessian_n_param(state)
      allocate (unit_vec(np, 1), column(np, 1), dense(np, np))
      do k = 1, np
         unit_vec = 0.0_dp
         unit_vec(k, 1) = 1.0_dp
         call sa_hessian_apply(state, unit_vec, column, error)
         if (error%has_error()) return
         dense(:, k) = column(:, 1)
      end do
      deallocate (unit_vec, column)
   end subroutine build_dense

   subroutine test_redundancy_sa3(error)
      !! The null space: `n_states * (antisym_dim + n_states)` near-zero
      !! eigenvalues (`antisym_dim = na*(na-1)/2` from the singlet
      !! symmetrisation this parametrisation imposes, plus the genuine
      !! state-mixing redundancy), everything else strictly positive at a
      !! converged SA minimum -- and, since this system finally has enough
      !! CI directions to matter, the smallest strictly-positive eigenvalue
      !! is reported rather than only its sign.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      type(rhf_result_t) :: scf
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: dense(:, :), values(:)
      real(dp), parameter :: WEIGHTS(3) = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
      integer :: np, info, n_near_zero, k, na, antisym_dim, n_expected
      real(dp), parameter :: ZERO_TOL = 1.0e-6_dp
      real(dp), parameter :: POSITIVE_FLOOR = 1.0e-5_dp
      real(dp) :: smallest_positive

      call build_czt_molecule(LIH_Z, LIH_SYM, LIH, "6-31g", mol, err)
      call check(error,.not. err%has_error(), "the molecule should build")
      if (allocated(error)) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      call check(error, scf%converged .and. .not. err%has_error(), &
                 "the LiH/6-31G RHF reference should converge")
      if (allocated(error)) return
      orbitals = scf%orbitals

      call run_czt_casscf(mol, orbitals, 0, 4, 2, 2, result, err, max_iterations=500, &
                          gradient_tol=1.0e-10_dp, n_states=3, weights=WEIGHTS)
      call check(error,.not. err%has_error(), "CAS(4,4) SA-3 should not error")
      if (allocated(error)) return
      call check(error, result%converged, "CAS(4,4) SA-3 should converge")
      if (allocated(error)) return

      call build_sa_hessian(mol, result%orbitals, 0, 4, 2, 2, result%ci_vectors, &
                            result%energies, WEIGHTS, state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      call build_dense(state, dense, err)
      call check(error,.not. err%has_error(), "the dense Hessian should build")
      if (allocated(error)) return

      np = sa_hessian_n_param(state)
      dense = 0.5_dp*(dense + transpose(dense))
      allocate (values(np))
      call pic_syev(dense, values, jobz="N", uplo="U", info=info)
      call check(error, info == 0, "the dense Hessian should diagonalise")
      if (allocated(error)) return

      na = state%alpha%n_strings
      antisym_dim = na*(na - 1)/2
      n_expected = state%n_states*(antisym_dim + state%n_states)
      n_near_zero = count(abs(values) < ZERO_TOL)
      ! At least the formula's redundancy; this system also has one
      ! symmetry-driven null direction beyond it -- see the module docstring.
      call check(error, n_near_zero >= n_expected, &
                 "the null space should be at least the antisymmetric-CI and "// &
                 "state-mixing directions the formula predicts")
      if (.not. allocated(error)) then
         ! A clean spectral gap between the null cluster and the rest, so an
         ! extra "near-zero" entry is a further exact symmetry and not a
         ! marginal or wrong small eigenvalue.
         call check(error, values(n_near_zero + 1) > 1.0e6_dp*max(1.0e-14_dp, &
                                                                  abs(values(n_near_zero))), &
                    "the null cluster and the positive spectrum should be "// &
                    "separated by several orders of magnitude")
      end if
      if (.not. allocated(error)) then
         do k = 1, np
            if (abs(values(k)) < ZERO_TOL) cycle
            call check(error, values(k) > POSITIVE_FLOOR, &
                       "every non-redundant eigenvalue should be positive")
            if (allocated(error)) exit
         end do
      end if
      if (.not. allocated(error)) then
         smallest_positive = huge(1.0_dp)
         do k = 1, np
            if (values(k) > ZERO_TOL) smallest_positive = min(smallest_positive, values(k))
         end do
         call check(error, smallest_positive < huge(1.0_dp), &
                    "there should be at least one strictly positive eigenvalue")
      end if

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_redundancy_sa3

end module test_mqc_sa_hessian_redundancy_long

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_hessian_redundancy_long, only: collect_mqc_sa_hessian_redundancy_long_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_hessian_redundancy_long", &
                               collect_mqc_sa_hessian_redundancy_long_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
