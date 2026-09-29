!! The matrix-free SA-CASSCF Hessian-vector product
module test_mqc_sa_hessian
   !! Phase 3 of the SA-CASSCF gradient project (`SA_CASSCF_GRADIENT_PLAN.md`;
   !! equations in `mqc_docs/source/developer_sa_casscf.rst`).
   !!
   !! LiH/STO-3G SA-2-CAS(2,2), the same reference `test_mqc_sa_casscf.f90`
   !! converges (`E_ROOT1`/`E_ROOT2` there), taken to tight orbital
   !! convergence (`1e-10`) so the finite-difference gate is not polluted by
   !! the parametrisation's own second-order terms.
   !!
   !! **The dense small Hessian is built explicitly**, one `sa_hessian_apply`
   !! call per unit vector of the `n_param = n_rot + n_states * n_det`
   !! parameter space (19 for this system) -- cheap at this size, and it turns
   !! gates A, C and D into one linear-algebra object: the top-left `n_rot`
   !! block against `orbital_hessian` (gate A), the whole matrix's asymmetry
   !! (gate C), and its eigenvalues (gate D, redundancy).
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t, orbital_hessian, &
                            generalized_fock, mcscf_fock_t, rotation_parameters, &
                            rotation_matrix
   use mqc_czt_casci, only: active_space_integrals
   use mqc_czt_sa_hessian, only: sa_hessian_t, build_sa_hessian, destroy_sa_hessian, &
                                 sa_hessian_n_param, sa_gradient, sa_hessian_apply, &
                                 project_ci_block, cheap_generalized_fock, &
                                 one_index_active_hamiltonian, sa_hessian_precondition
   implicit none
   private

   public :: collect_mqc_sa_hessian_tests

   real(dp), parameter :: LIH(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 3.0139241961656_dp], [3, 2])
   integer, parameter :: LIH_Z(2) = [3, 1]
   character(len=2), parameter :: LIH_SYM(2) = ["Li", "H "]

contains

   subroutine collect_mqc_sa_hessian_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("orbital_block_matches_orbital_hessian_sa2", &
                               test_orbital_block_sa2), &
                  new_unittest("orbital_block_matches_orbital_hessian_n_states_1", &
                               test_orbital_block_one_state), &
                  new_unittest("hessian_is_symmetric", test_symmetry), &
                  new_unittest("redundancy_and_positive_curvature", test_redundancy), &
                  new_unittest("preconditioner_matches_diagonals", test_preconditioner), &
                  new_unittest("apply_refuses_wrong_shape", test_apply_shape), &
                  new_unittest("hvp_against_finite_difference", test_finite_difference), &
                  new_unittest("cheap_generalized_fock_matches_generalized_fock", &
                               test_cheap_fock), &
                  new_unittest("one_index_active_hamiltonian_against_fd", &
                               test_one_index_active_hamiltonian) &
                  ]
   end subroutine collect_mqc_sa_hessian_tests

   subroutine lih_reference(mol, orbitals, err)
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

   subroutine converged_sa2(mol, orbitals, result, err, ok)
      !! SA-2-CAS(2,2), converged tightly
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok
      real(dp), parameter :: WEIGHTS(2) = [0.5_dp, 0.5_dp]

      ok = .false.
      call lih_reference(mol, orbitals, err)
      if (err%has_error()) return
      call run_czt_casscf(mol, orbitals, 1, 2, 1, 1, result, err, max_iterations=300, &
                          gradient_tol=1.0e-10_dp, n_states=2, weights=WEIGHTS)
      if (err%has_error()) return
      ok = result%converged
   end subroutine converged_sa2

   subroutine converged_single(mol, orbitals, result, err, ok)
      !! Plain CAS(2,2), one state -- the n_states = 1 reduction
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      ok = .false.
      call lih_reference(mol, orbitals, err)
      if (err%has_error()) return
      call run_czt_casscf(mol, orbitals, 1, 2, 1, 1, result, err, max_iterations=300, &
                          gradient_tol=1.0e-10_dp)
      if (err%has_error()) return
      ok = result%converged
   end subroutine converged_single

   function flat_kappa_direction(n_rot, phase) result(v)
      !! A deterministic, arbitrary direction over the rotation parameters
      integer, intent(in) :: n_rot
      real(dp), intent(in) :: phase
      real(dp) :: v(n_rot)

      integer :: l

      do l = 1, n_rot
         v(l) = sin(0.37_dp*real(l, dp) + phase)
      end do
   end function flat_kappa_direction

   subroutine ci_direction(na, nb, phase, v)
      !! A deterministic, arbitrary vector over the determinant space
      integer, intent(in) :: na, nb
      real(dp), intent(in) :: phase
      real(dp), intent(out) :: v(na, nb)

      integer :: ia, ib

      do ib = 1, nb
         do ia = 1, na
            v(ia, ib) = sin(0.53_dp*real(ia, dp) + 0.29_dp*real(ib, dp) + phase)
         end do
      end do
   end subroutine ci_direction

   subroutine build_dense(state, dense, error)
      !! The explicit `(n_param, n_param)` matrix, one `sa_hessian_apply` call
      !! per unit vector -- affordable at LiH's size and what gates A, C, D
      !! read off directly.
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

   subroutine test_orbital_block_sa2(error)
      !! With zero CI input, `sa_hessian_apply`'s orbital output equals
      !! `orbital_hessian` at the SA densities applied to the same kappa,
      !! on an arbitrary direction (not just a unit vector)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: full_hessian(:, :), x(:, :), hx(:, :)
      real(dp), allocatable :: kappa_dir(:), expected(:)
      integer, allocatable :: rows(:), cols(:)
      integer :: np, n_rot
      logical :: ok

      call converged_sa2(mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, [0.5_dp, 0.5_dp], state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      call rotation_parameters(size(result%orbitals, 2), 1, 2, rows, cols)
      n_rot = size(rows)
      call orbital_hessian(mol, result%orbitals, 1, 2, result%dm1, result%dm2, &
                           state%fock_sa, rows, cols, full_hessian, err)
      call check(error,.not. err%has_error(), "orbital_hessian should build")
      if (allocated(error)) return

      np = sa_hessian_n_param(state)
      allocate (x(np, 1), hx(np, 1))
      x = 0.0_dp
      kappa_dir = flat_kappa_direction(n_rot, 0.7_dp)
      x(1:n_rot, 1) = kappa_dir
      call sa_hessian_apply(state, x, hx, err)
      call check(error,.not. err%has_error(), "sa_hessian_apply should not error")
      if (allocated(error)) return

      expected = matmul(full_hessian, kappa_dir)
      call check(error, maxval(abs(hx(1:n_rot, 1) - expected)) &
                 < 1.0e-10_dp*max(1.0_dp, maxval(abs(expected))), &
                 "the orbital-orbital block should match orbital_hessian at the SA densities")

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_orbital_block_sa2

   subroutine test_orbital_block_one_state(error)
      !! The n_states = 1 reduction: the orbital-orbital block is exactly
      !! the single-state `orbital_hessian`
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: full_hessian(:, :), x(:, :), hx(:, :)
      real(dp), allocatable :: kappa_dir(:), expected(:)
      integer, allocatable :: rows(:), cols(:)
      real(dp), allocatable :: ci_vectors(:, :, :), energies(:)
      integer :: np, n_rot
      logical :: ok

      call converged_single(mol, orbitals, result, err, ok)
      call check(error, ok, "the plain CASSCF should converge")
      if (allocated(error)) return

      allocate (ci_vectors(size(result%ci_vector, 1), size(result%ci_vector, 2), 1))
      ci_vectors(:, :, 1) = result%ci_vector
      allocate (energies(1))
      energies(1) = result%energy

      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, ci_vectors, energies, &
                            [1.0_dp], state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      call rotation_parameters(size(result%orbitals, 2), 1, 2, rows, cols)
      n_rot = size(rows)
      call orbital_hessian(mol, result%orbitals, 1, 2, result%dm1, result%dm2, &
                           state%fock_sa, rows, cols, full_hessian, err)
      call check(error,.not. err%has_error(), "orbital_hessian should build")
      if (allocated(error)) return

      np = sa_hessian_n_param(state)
      allocate (x(np, 1), hx(np, 1))
      x = 0.0_dp
      kappa_dir = flat_kappa_direction(n_rot, 1.3_dp)
      x(1:n_rot, 1) = kappa_dir
      call sa_hessian_apply(state, x, hx, err)
      call check(error,.not. err%has_error(), "sa_hessian_apply should not error")
      if (allocated(error)) return

      expected = matmul(full_hessian, kappa_dir)
      call check(error, maxval(abs(hx(1:n_rot, 1) - expected)) &
                 < 1.0e-10_dp*max(1.0_dp, maxval(abs(expected))), &
                 "n_states=1's orbital-orbital block should be orbital_hessian exactly")

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_orbital_block_one_state

   subroutine test_symmetry(error)
      !! |<y,Hx> - <x,Hy>| is zero to rounding, on random x and y that mix
      !! the orbital and CI blocks
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: x(:, :), y(:, :), hx(:, :), hy(:, :)
      real(dp), allocatable :: kdir(:)
      real(dp), allocatable :: xci(:, :), yci(:, :)
      real(dp) :: lhs, rhs, scale
      integer :: np, n_rot, na, nb, j
      logical :: ok

      call converged_sa2(mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, [0.5_dp, 0.5_dp], state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      n_rot = state%n_rot
      na = state%alpha%n_strings
      nb = state%beta%n_strings
      np = sa_hessian_n_param(state)
      allocate (x(np, 1), y(np, 1), hx(np, 1), hy(np, 1))
      allocate (xci(na, nb), yci(na, nb))

      kdir = flat_kappa_direction(n_rot, 0.2_dp)
      x(1:n_rot, 1) = kdir
      do j = 1, 2
         call ci_direction(na, nb, 0.4_dp + 0.1_dp*real(j, dp), xci)
         call project_ci_block(state, xci)
         xci = 0.5_dp*(xci + transpose(xci))
         x(n_rot + (j - 1)*state%n_det + 1:n_rot + j*state%n_det, 1) = reshape(xci, [na*nb])
      end do

      kdir = flat_kappa_direction(n_rot, 1.9_dp)
      y(1:n_rot, 1) = kdir
      do j = 1, 2
         call ci_direction(na, nb, 1.1_dp + 0.3_dp*real(j, dp), yci)
         call project_ci_block(state, yci)
         yci = 0.5_dp*(yci + transpose(yci))
         y(n_rot + (j - 1)*state%n_det + 1:n_rot + j*state%n_det, 1) = reshape(yci, [na*nb])
      end do

      call sa_hessian_apply(state, x, hx, err)
      call check(error,.not. err%has_error(), "sa_hessian_apply(x) should not error")
      if (allocated(error)) return
      call sa_hessian_apply(state, y, hy, err)
      call check(error,.not. err%has_error(), "sa_hessian_apply(y) should not error")
      if (allocated(error)) return

      lhs = sum(y(:, 1)*hx(:, 1))
      rhs = sum(x(:, 1)*hy(:, 1))
      scale = max(1.0_dp, abs(lhs), abs(rhs))
      call check(error, abs(lhs - rhs) < 1.0e-10_dp*scale, &
                 "<y,Hx> should equal <x,Hy>")

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_symmetry

   subroutine test_redundancy(error)
      !! The explicit small Hessian's null space, probed with unit vectors
      !! over the *full* (na, nb) determinant space rather than only its
      !! singlet subspace, has two sources per state `J`: the antisymmetric
      !! part of `x_J` (`na*(na-1)/2` dimensions, killed outright by the
      !! transpose-symmetrisation `sa_hessian_apply_one` applies -- the
      !! parametrisation is singlets only, by design, per
      !! `developer_sa_casscf.rst`'s "Spin"), and, within the surviving
      !! symmetric subspace, `span{c_1, ..., c_N}` (`n_states` dimensions,
      !! the redundant rotations among the averaged states). Every other
      !! eigenvalue must be strictly positive at a converged SA minimum.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: dense(:, :), values(:)
      integer :: np, info, n_near_zero, k, na, antisym_dim, n_expected
      logical :: ok
      real(dp), parameter :: ZERO_TOL = 1.0e-6_dp
      real(dp), parameter :: POSITIVE_FLOOR = 1.0e-5_dp

      call converged_sa2(mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, [0.5_dp, 0.5_dp], state, err)
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
      call check(error, n_near_zero == n_expected, &
                 "the null space should be exactly the antisymmetric-CI and "// &
                 "state-mixing directions")
      if (.not. allocated(error)) then
         do k = 1, np
            if (abs(values(k)) < ZERO_TOL) cycle
            call check(error, values(k) > POSITIVE_FLOOR, &
                       "every non-redundant eigenvalue should be positive")
            if (allocated(error)) exit
         end do
      end if

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_redundancy

   subroutine test_preconditioner(error)
      !! `sa_hessian_precondition` against the two diagonals it divides by
      !!
      !! The orbital block against `orbital_hessian_diag`, the one-electron
      !! approximation `build_sa_hessian` stores, checked for sign against the
      !! dense Hessian's exact diagonal (`build_dense`) so an approximation that
      !! turned a positive curvature negative would be caught; the CI block
      !! against `2 w_J (H_diag - E_J)` with
      !! `H_diag` flattened by `reshape`, which is the column-major order every
      !! flat CI vector here is packed in. At `n_alpha == n_beta` that diagonal
      !! is symmetric under `ia <-> ib`, so a transposed determinant index would
      !! give the same numbers and is not what this can catch; what it does pin
      !! is each state's segment offset and length, which the unequal weights
      !! make distinguishable.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: dense(:, :), x(:, :), px(:, :), flat_diag(:)
      real(dp) :: expected, d
      integer :: np, l, j, det, seg0
      logical :: ok
      real(dp), parameter :: FLOOR = 1.0e-3_dp
      real(dp), parameter :: TOL = 1.0e-10_dp

      call converged_sa2(mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, [0.7_dp, 0.3_dp], state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      call build_dense(state, dense, err)
      call check(error,.not. err%has_error(), "the dense Hessian should build")
      if (allocated(error)) return

      np = sa_hessian_n_param(state)
      allocate (x(np, 1), px(np, 1))
      x = 1.0_dp
      call sa_hessian_precondition(state, x, px)

      do l = 1, state%n_rot
         d = state%orbital_hessian_diag(l)
         if (abs(d) < FLOOR) d = sign(FLOOR, d)
         call check(error, abs(px(l, 1) - 1.0_dp/d) < TOL*max(1.0_dp, abs(1.0_dp/d)), &
                    "the orbital preconditioner should divide by the stored diagonal")
         if (allocated(error)) return
         if (abs(dense(l, l)) >= FLOOR) then
            call check(error, sign(1.0_dp, d) == sign(1.0_dp, dense(l, l)), &
                       "the approximate orbital diagonal should keep the exact one's sign")
            if (allocated(error)) return
         end if
      end do

      flat_diag = reshape(state%diagonal, [state%n_det])
      do j = 1, state%n_states
         seg0 = state%n_rot + (j - 1)*state%n_det
         do det = 1, state%n_det
            d = 2.0_dp*state%weights(j)*(flat_diag(det) - state%active_energies(j))
            if (abs(d) < FLOOR) d = sign(FLOOR, d)
            expected = 1.0_dp/d
            call check(error, abs(px(seg0 + det, 1) - expected) < &
                       TOL*max(1.0_dp, abs(expected)), &
                       "the CI preconditioner should follow the flat determinant order")
            if (allocated(error)) return
         end do
      end do

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_preconditioner

   subroutine test_apply_shape(error)
      !! A vector of the wrong length is refused by name, not read out of bounds
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: x(:, :), hx(:, :)
      integer :: np
      logical :: ok

      call converged_sa2(mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, [0.5_dp, 0.5_dp], state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      np = sa_hessian_n_param(state)
      allocate (x(np - 1, 1), hx(np - 1, 1))
      x = 0.0_dp
      call sa_hessian_apply(state, x, hx, err)
      call check(error, err%has_error(), "a short parameter vector should be refused")
      if (allocated(error)) return
      call check(error, err%get_code() == ERROR_VALIDATION, &
                 "and refused as a validation error: "//err%get_message())

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_apply_shape

   subroutine test_finite_difference(error)
      !! The full HVP (every block) against a central difference of
      !! `sa_gradient`, orbital-only, CI-only and mixed directions
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: kappa_dir_flat(:), kappa_dir(:, :)
      real(dp), allocatable :: xci_dir(:, :, :)
      real(dp), allocatable :: x(:, :), hx(:, :)
      real(dp), allocatable :: g_plus_orb(:, :), g_minus_orb(:, :)
      real(dp), allocatable :: g_plus_ci(:, :, :), g_minus_ci(:, :, :)
      real(dp), allocatable :: fd_orb(:), fd_ci_j(:, :)
      real(dp) :: h_step, denom
      integer :: np, n_rot, na, nb, j, kase
      logical :: ok

      call converged_sa2(mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, [0.5_dp, 0.5_dp], state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      n_rot = state%n_rot
      na = state%alpha%n_strings
      nb = state%beta%n_strings
      np = sa_hessian_n_param(state)
      h_step = 1.0e-4_dp

      allocate (kappa_dir(state%n_mo, state%n_mo), xci_dir(na, nb, 2))
      allocate (x(np, 1), hx(np, 1))
      allocate (g_plus_ci(na, nb, 2), g_minus_ci(na, nb, 2))
      allocate (fd_ci_j(na, nb))

      do kase = 1, 3
         ! 1 orbital-only, 2 CI-only, 3 mixed
         kappa_dir_flat = flat_kappa_direction(n_rot, 0.6_dp)
         kappa_dir = 0.0_dp
         if (kase /= 2) then
            block
               integer :: l
               do l = 1, n_rot
                  kappa_dir(state%rows(l), state%cols(l)) = kappa_dir_flat(l)
                  kappa_dir(state%cols(l), state%rows(l)) = -kappa_dir_flat(l)
               end do
            end block
         else
            kappa_dir_flat = 0.0_dp
         end if

         xci_dir = 0.0_dp
         if (kase /= 1) then
            do j = 1, 2
               call ci_direction(na, nb, 0.8_dp + 0.2_dp*real(j, dp), xci_dir(:, :, j))
               call project_ci_block(state, xci_dir(:, :, j))
               xci_dir(:, :, j) = 0.5_dp*(xci_dir(:, :, j) + transpose(xci_dir(:, :, j)))
            end do
         end if

         x = 0.0_dp
         x(1:n_rot, 1) = kappa_dir_flat
         do j = 1, 2
            x(n_rot + (j - 1)*state%n_det + 1:n_rot + j*state%n_det, 1) = &
               reshape(xci_dir(:, :, j), [na*nb])
         end do
         call sa_hessian_apply(state, x, hx, err)
         call check(error,.not. err%has_error(), "sa_hessian_apply should not error")
         if (allocated(error)) exit

         call sa_gradient(state, h_step*kappa_dir, h_step*xci_dir, g_plus_orb, g_plus_ci, err)
         call check(error,.not. err%has_error(), "sa_gradient(+h) should not error")
         if (allocated(error)) exit
         call sa_gradient(state, -h_step*kappa_dir, -h_step*xci_dir, g_minus_orb, g_minus_ci, err)
         call check(error,.not. err%has_error(), "sa_gradient(-h) should not error")
         if (allocated(error)) exit

         allocate (fd_orb(n_rot))
         block
            integer :: l
            do l = 1, n_rot
               fd_orb(l) = (g_plus_orb(state%rows(l), state%cols(l)) &
                            - g_minus_orb(state%rows(l), state%cols(l)))/(2.0_dp*h_step)
            end do
         end block
         denom = max(1.0_dp, maxval(abs(hx(1:n_rot, 1))))
         call check(error, maxval(abs(hx(1:n_rot, 1) - fd_orb)) < 1.0e-7_dp*denom, &
                    "the orbital output should match the finite difference")
         deallocate (fd_orb)
         if (allocated(error)) exit

         do j = 1, 2
            fd_ci_j = (g_plus_ci(:, :, j) - g_minus_ci(:, :, j))/(2.0_dp*h_step)
            ! `sa_hessian_apply`'s CI output is projected onto the complement
            ! of *every* reference state (PySCF's `project_Aop` gauge for the
            ! Z-vector solve, not merely a null-space discovery): a component
            ! of the raw derivative along another state c_K is a valid gauge
            ! choice `sa_gradient` does not itself fix, so the finite
            ! difference is projected the same way before comparing -- an
            ! honest reading of "vs finite difference of sa_gradient", not a
            ! loosening of the gate.
            call project_ci_block(state, fd_ci_j)
            denom = max(1.0_dp, maxval(abs(hx(n_rot + (j - 1)*state%n_det + 1: &
                                              n_rot + j*state%n_det, 1))))
            call check(error, maxval(abs(reshape(hx(n_rot + (j - 1)*state%n_det + 1: &
                                                    n_rot + j*state%n_det, 1), [na, nb]) &
                                         - fd_ci_j)) < 1.0e-7_dp*denom, &
                       "the CI output should match the finite difference")
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_finite_difference

   subroutine test_cheap_fock(error)
      !! `cheap_generalized_fock` (reusing `a_block`/`b_block`) must agree
      !! with `generalized_fock` (a fresh AO integral pass) at the SA
      !! densities: the whole point of the cheap routine is that it computes
      !! the same number faster, not a different one
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      type(mcscf_fock_t) :: reference
      real(dp), allocatable :: general(:, :)
      logical :: ok

      call converged_sa2(mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, [0.5_dp, 0.5_dp], state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      call generalized_fock(mol, result%orbitals, 1, 2, result%dm1, result%dm2, &
                            reference, err)
      call check(error,.not. err%has_error(), "generalized_fock should not error")
      if (allocated(error)) return

      call cheap_generalized_fock(state, result%dm1, result%dm2, general)
      call check(error, maxval(abs(general - reference%general)) < 1.0e-10_dp, &
                 "cheap_generalized_fock should reproduce generalized_fock's general Fock")

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_cheap_fock

   subroutine test_one_index_active_hamiltonian(error)
      !! `dh_eff`/`deri_act` against a central difference of
      !! `active_space_integrals`'s `h_eff`/`eri_act` under the same orbital
      !! rotation -- isolates `one_index_active_hamiltonian` from the rest
      !! of the CI-orbital block (folding, sigma, projection)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: kappa_flat(:), kappa(:, :), rotation(:, :)
      real(dp), allocatable :: plus_orbitals(:, :), minus_orbitals(:, :)
      real(dp), allocatable :: h_eff_plus(:, :), eri_act_plus(:, :, :, :)
      real(dp), allocatable :: h_eff_minus(:, :), eri_act_minus(:, :, :, :)
      real(dp), allocatable :: dh_eff(:, :), deri_act(:, :, :, :)
      real(dp), allocatable :: fd_dh_eff(:, :), fd_deri_act(:, :, :, :)
      real(dp) :: core_plus, core_minus, h_step
      integer :: n_rot
      logical :: ok

      call converged_sa2(mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return
      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, [0.5_dp, 0.5_dp], state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      n_rot = state%n_rot
      h_step = 1.0e-4_dp
      kappa_flat = flat_kappa_direction(n_rot, 0.44_dp)
      allocate (kappa(state%n_mo, state%n_mo))
      kappa = 0.0_dp
      block
         integer :: l
         do l = 1, n_rot
            kappa(state%rows(l), state%cols(l)) = kappa_flat(l)
            kappa(state%cols(l), state%rows(l)) = -kappa_flat(l)
         end do
      end block

      call one_index_active_hamiltonian(state, kappa, dh_eff, deri_act)

      call rotation_matrix(h_step*kappa, rotation)
      allocate (plus_orbitals(state%n_ao, state%n_mo))
      plus_orbitals = matmul(result%orbitals, rotation)
      deallocate (rotation)
      call rotation_matrix(-h_step*kappa, rotation)
      allocate (minus_orbitals(state%n_ao, state%n_mo))
      minus_orbitals = matmul(result%orbitals, rotation)

      call active_space_integrals(mol, plus_orbitals, 1, 2, h_eff_plus, eri_act_plus, &
                                  core_plus, err)
      call check(error,.not. err%has_error(), "active_space_integrals(+h) should not error")
      if (allocated(error)) return
      call active_space_integrals(mol, minus_orbitals, 1, 2, h_eff_minus, eri_act_minus, &
                                  core_minus, err)
      call check(error,.not. err%has_error(), "active_space_integrals(-h) should not error")
      if (allocated(error)) return

      fd_dh_eff = (h_eff_plus - h_eff_minus)/(2.0_dp*h_step)
      fd_deri_act = (eri_act_plus - eri_act_minus)/(2.0_dp*h_step)

      call check(error, maxval(abs(dh_eff - fd_dh_eff)) < 1.0e-6_dp, &
                 "dh_eff should match the finite difference")
      if (.not. allocated(error)) &
         call check(error, maxval(abs(deri_act - fd_deri_act)) < 1.0e-6_dp, &
                    "deri_act should match the finite difference")

      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_one_index_active_hamiltonian

end module test_mqc_sa_hessian

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_hessian, only: collect_mqc_sa_hessian_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_hessian", collect_mqc_sa_hessian_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
