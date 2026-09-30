!! The single-root SA-CASSCF nuclear gradient via the Z-vector
module test_mqc_sa_gradient
   !! Phase 4 of the SA-CASSCF gradient project (`SA_CASSCF_GRADIENT_PLAN.md`;
   !! equations in `mqc_docs/source/developer_sa_casscf.rst`).
   !!
   !! LiH/STO-3G SA-2-CAS(2,2) is `test_mqc_sa_hessian.f90`'s system, reused
   !! here for the same reason: fast, and every gate that does not need a
   !! genuinely multi-dimensional CI space fits on it. LiH/6-31G CAS(4,4)
   !! SA-3 (`test_mqc_sa_hessian_multistate.f90`'s system) exercises the
   !! CI-response and orbital-CI pieces with 18 non-redundant CI directions
   !! per state rather than one.
   !!
   !! **The reference used for the per-piece gates matches the connection the
   !! analytic code assumes.** `orbital_response_gradient`/`ci_response_gradient`
   !! both contain an overlap (Pulay) term, valid for the same reason
   !! `czt_mcscf_gradient`'s own is: as R moves, a *fixed numerical* MO
   !! coefficient matrix `C0` is re-expressed in the symmetric-orthogonalised
   !! connection `C(R) = C0 (C0^T S(R) C0)^(-1/2)`, not held as `C0` verbatim
   !! against a now-non-orthonormal metric. `connection_orbitals` builds
   !! exactly this, and `root_energy_at`/`orbital_response_phi_at`/
   !! `ci_response_phi_at` central-difference the three pieces of the
   !! Lagrangian (`E_I`, `kappa_bar . dE_SA/dkappa`, `sum_J xbar_J . dE_J/dc_J`)
   !! at fixed `kappa_bar`/`xbar`/CI vectors in that connection -- an
   !! independent reference an earlier version of this file did not have: it
   !! held `C(R) = C0` verbatim, a *different* connection than the analytic
   !! code's overlap term assumes, and so appeared to show that term was
   !! unneeded.
   !!
   !! **Gate 2 (the weighted-sum identity) is structurally weak and cannot
   !! replace gate 3.** Because the Z-vector equation is linear in its
   !! right-hand side and `sum_I w_I RHS_I = RHS(D_SA, d_SA) = 0` at
   !! convergence, `sum_I w_I kappa_bar_I = 0` and `sum_I w_I xbar_I = 0`
   !! *exactly*, by linearity alone -- for any number of states, regardless
   !! of whether the response formula contracting them is correct, as long as
   !! it is linear/homogeneous in `(kappa_bar, xbar)` (which any plausible
   !! response formula is). A response piece that is scaled or mis-signed
   !! uniformly still passes gate 2. It is kept below because it is cheap and
   !! does check the base term and the RHS/solve plumbing, but every response
   !! formula change in this project has needed gate 3 (or the per-piece
   !! connection references below) to actually catch anything.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use omp_lib, only: omp_get_max_threads, omp_set_num_threads
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_lapack_interfaces, only: pic_syev
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, build_czt_molecule
   use mqc_czt_rhf, only: rhf_result_t, run_czt_rhf
   use mqc_czt_mcscf, only: run_czt_casscf, casscf_result_t, sa_density_matrices, &
                            generalized_fock, orbital_gradient, mcscf_fock_t, &
                            rotation_parameters
   use mqc_czt_mcscf_gradient, only: czt_mcscf_gradient
   use mqc_czt_sa_gradient, only: czt_sa_casscf_gradient, sa_casscf_gradient_general, &
                                  orbital_response_gradient, ci_response_gradient, &
                                  sa_zvector_solve
   use mqc_czt_sa_hessian, only: sa_hessian_t, build_sa_hessian, destroy_sa_hessian, &
                                 sa_hessian_n_param, project_ci_block
   use mqc_czt_casci, only: active_space_integrals
   use mqc_ci, only: absorb_one_electron, sigma_vector
   use mqc_rdm, only: active_space_rdms, rdm_energy
   use mqc_determinants, only: link_table_t, build_link_table
   implicit none
   private

   public :: collect_mqc_sa_gradient_tests

   real(dp), parameter :: LIH(3, 2) = reshape( &
                          [0.0_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 3.0139241961656_dp], [3, 2])
   integer, parameter :: LIH_Z(2) = [3, 1]
   character(len=2), parameter :: LIH_SYM(2) = ["Li", "H "]
   real(dp), parameter :: SA2_WEIGHTS(2) = [0.5_dp, 0.5_dp]

contains

   subroutine collect_mqc_sa_gradient_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("n_states_1_is_bit_identical_to_czt_mcscf_gradient", &
                               test_short_circuit_bit_identical), &
                  new_unittest("general_zvector_path_at_n_states_1_agrees", &
                               test_general_path_n1), &
                  new_unittest("weighted_sum_of_root_gradients_is_the_sa_gradient", &
                               test_weighted_sum_identity), &
                  new_unittest("each_root_gradient_sums_to_zero_over_atoms", &
                               test_translational_invariance), &
                  new_unittest("unequal_weights_are_refused", test_unequal_weights_refused), &
                  new_unittest("weighted_sum_holds_on_cas44_sa3", &
                               test_weighted_sum_cas44_sa3), &
                  new_unittest("base_term_matches_connection_reference", &
                               test_base_term_connection), &
                  new_unittest("orbital_response_matches_connection_reference", &
                               test_orbital_response_connection), &
                  new_unittest("ci_response_matches_connection_reference", &
                               test_ci_response_connection), &
                  new_unittest("root_gradients_difference_their_own_energies", &
                               test_finite_difference_lih_sa2), &
                  new_unittest("cas44_sa3_root_gradients_difference_their_energies", &
                               test_finite_difference_cas44_sa3) &
                  ]
   end subroutine collect_mqc_sa_gradient_tests

   ! ============================================================
   ! Setup helpers
   ! ============================================================

   subroutine lih_reference(mol, orbitals, basis, err)
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      character(len=*), intent(in) :: basis
      type(error_t), intent(inout) :: err

      type(rhf_result_t) :: scf

      call build_czt_molecule(LIH_Z, LIH_SYM, LIH, basis, mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error()) return
      if (.not. scf%converged) then
         call err%set(ERROR_VALIDATION, "the LiH RHF reference did not converge")
         return
      end if
      orbitals = scf%orbitals
   end subroutine lih_reference

   subroutine converged_sa(coordinates, basis, n_inactive, n_active, weights, mol, &
                           orbitals, result, err, ok)
      !! A converged SA-CASSCF at an arbitrary geometry, `size(weights)` states
      real(dp), intent(in) :: coordinates(3, 2)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: weights(:)
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok

      integer :: n_half

      type(rhf_result_t) :: scf

      ok = .false.
      call build_czt_molecule(LIH_Z, LIH_SYM, coordinates, basis, mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error()) return
      if (.not. scf%converged) then
         call err%set(ERROR_VALIDATION, "the LiH RHF reference did not converge")
         return
      end if
      ! LiH has four electrons; whatever the inactive space leaves is active.
      n_half = (4 - 2*n_inactive)/2
      if (size(weights) == 1) then
         call run_czt_casscf(mol, scf%orbitals, n_inactive, n_active, n_half, n_half, result, err, &
                             max_iterations=400, gradient_tol=1.0e-10_dp)
      else
         call run_czt_casscf(mol, scf%orbitals, n_inactive, n_active, n_half, n_half, result, err, &
                             max_iterations=400, gradient_tol=1.0e-10_dp, &
                             n_states=size(weights), weights=weights)
      end if
      if (err%has_error()) return
      orbitals = result%orbitals
      ok = result%converged
   end subroutine converged_sa

   subroutine connection_orbitals(c0, mol, c_r, err)
      !! `C(R) = C0 (C0^T S(R) C0)^(-1/2)`: the symmetric-orthogonalisation
      !! connection that keeps a *fixed* numerical coefficient matrix `C0`
      !! re-orthonormalised against the AO metric `S(R)` at a displaced
      !! geometry, without touching the numbers in `C0` itself beyond that
      !! renormalisation. `czt_mcscf_gradient`'s own overlap (Pulay) term is
      !! the derivative of exactly this connection at `R0`, so it -- not
      !! `C(R) = C0` verbatim -- is the reference frame any Lagrange term
      !! built the same way needs to be checked in.
      real(dp), intent(in) :: c0(:, :)         !! (n_ao, n_mo), orthonormal at the reference R
      type(czt_molecule_t), intent(in) :: mol   !! At the (possibly displaced) geometry
      real(dp), allocatable, intent(out) :: c_r(:, :)
      type(error_t), intent(inout) :: err

      real(dp), allocatable :: s(:, :), m(:, :), work(:, :), minvsqrt(:, :)
      real(dp), allocatable :: eigvals(:)
      integer :: n_ao, n_mo, i, info

      if (err%has_error()) return
      n_ao = size(c0, 1)
      n_mo = size(c0, 2)

      call mol%overlap(s)
      allocate (work(n_ao, n_mo))
      call pic_gemm(s, c0, work, beta=0.0_dp)
      allocate (m(n_mo, n_mo))
      call pic_gemm(c0, work, m, transa="T", beta=0.0_dp)
      deallocate (work)

      allocate (eigvals(n_mo))
      call pic_syev(m, eigvals, jobz="V", uplo="U", info=info)
      if (info /= 0 .or. any(eigvals <= 0.0_dp)) then
         call err%set(ERROR_VALIDATION, "connection_orbitals: C0^T S(R) C0 is not "// &
                      "positive definite -- the displacement is too large")
         return
      end if

      allocate (work(n_mo, n_mo), minvsqrt(n_mo, n_mo))
      do i = 1, n_mo
         work(:, i) = m(:, i)/sqrt(eigvals(i))
      end do
      call pic_gemm(work, m, minvsqrt, transb="T", beta=0.0_dp)
      deallocate (work, m, eigvals)

      allocate (c_r(n_ao, n_mo))
      call pic_gemm(c0, minvsqrt, c_r, beta=0.0_dp)
   end subroutine connection_orbitals

   subroutine root_energy_at(mol_coords, basis, c0, n_inactive, n_active, dm1_i, dm2_i, &
                             energy, err)
      !! `E_I(C(R), c_I)` at fixed `c0`/`dm1_i`/`dm2_i`, in the connection
      real(dp), intent(in) :: mol_coords(3, 2)
      character(len=*), intent(in) :: basis
      real(dp), intent(in) :: c0(:, :)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1_i(:, :), dm2_i(:, :, :, :)
      real(dp), intent(out) :: energy
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      real(dp), allocatable :: c_r(:, :), h_eff(:, :), eri_act(:, :, :, :)
      real(dp) :: core_energy

      call build_czt_molecule(LIH_Z, LIH_SYM, mol_coords, basis, mol, err)
      if (err%has_error()) return
      call connection_orbitals(c0, mol, c_r, err)
      if (err%has_error()) return
      call active_space_integrals(mol, c_r, n_inactive, n_active, h_eff, eri_act, &
                                  core_energy, err)
      if (err%has_error()) return
      energy = core_energy + rdm_energy(h_eff, eri_act, dm1_i, dm2_i)
      call mol%destroy()
   end subroutine root_energy_at

   subroutine orbital_response_phi_at(mol_coords, basis, c0, n_inactive, n_active, &
                                      dm1_sa, dm2_sa, kappa_flat, rows, cols, phi, err)
      !! `kappa_bar . dE_SA/dkappa` at `C(R)`, fixed `dm1_sa`/`dm2_sa`
      real(dp), intent(in) :: mol_coords(3, 2)
      character(len=*), intent(in) :: basis
      real(dp), intent(in) :: c0(:, :)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1_sa(:, :), dm2_sa(:, :, :, :)
      real(dp), intent(in) :: kappa_flat(:)
      integer, intent(in) :: rows(:), cols(:)
      real(dp), intent(out) :: phi
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      real(dp), allocatable :: c_r(:, :), g_full(:, :)
      type(mcscf_fock_t) :: fock
      integer :: l

      call build_czt_molecule(LIH_Z, LIH_SYM, mol_coords, basis, mol, err)
      if (err%has_error()) return
      call connection_orbitals(c0, mol, c_r, err)
      if (err%has_error()) return
      call generalized_fock(mol, c_r, n_inactive, n_active, dm1_sa, dm2_sa, fock, err)
      if (err%has_error()) return
      call orbital_gradient(fock, n_inactive, n_active, g_full)
      phi = 0.0_dp
      do l = 1, size(rows)
         phi = phi + kappa_flat(l)*g_full(rows(l), cols(l))
      end do
      call mol%destroy()
   end subroutine orbital_response_phi_at

   subroutine ci_response_phi_at(mol_coords, basis, c0, n_inactive, n_active, n_alpha, &
                                 n_beta, ci_vectors, weights, xbar, alpha, beta, phi, err)
      !! `sum_J w_J xbar_J . dE_J/dc_J` at `C(R)`, fixed `ci_vectors`/`xbar`
      !!
      !! `dE_J/dc_J`, at `kappa=0` and evaluated *at* `c_J` (not `c_J+xbar_J`:
      !! `xbar_J` is a fixed Lagrange multiplier dotted with the gradient
      !! here, not a displacement of the point the gradient is taken at), is
      !! `2 w_J (H(C(R)) - E_J) c_J` with `E_J = <c_J|H(C(R))|c_J>` -- the
      !! same construction `sa_gradient` (`mqc_czt_sa_hessian.f90`) uses at
      !! its own `x = 0`.
      real(dp), intent(in) :: mol_coords(3, 2)
      character(len=*), intent(in) :: basis
      real(dp), intent(in) :: c0(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)     !! (na, nb, n_states)
      real(dp), intent(in) :: weights(:)
      real(dp), intent(in) :: xbar(:, :, :)           !! (na, nb, n_states)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), intent(out) :: phi
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      real(dp), allocatable :: c_r(:, :), h_eff(:, :), eri_act(:, :, :, :), folded(:, :)
      real(dp), allocatable :: sigma_j(:, :), gradient_ci_j(:, :)
      real(dp) :: core_energy, e_active
      integer :: j, n_states

      call build_czt_molecule(LIH_Z, LIH_SYM, mol_coords, basis, mol, err)
      if (err%has_error()) return
      call connection_orbitals(c0, mol, c_r, err)
      if (err%has_error()) return
      call active_space_integrals(mol, c_r, n_inactive, n_active, h_eff, eri_act, &
                                  core_energy, err)
      if (err%has_error()) return
      call absorb_one_electron(h_eff, eri_act, n_alpha + n_beta, folded, err)
      if (err%has_error()) return

      n_states = size(weights)
      allocate (sigma_j(alpha%n_strings, beta%n_strings))
      allocate (gradient_ci_j(alpha%n_strings, beta%n_strings))
      phi = 0.0_dp
      do j = 1, n_states
         call sigma_vector(folded, ci_vectors(:, :, j), alpha, beta, sigma_j, err)
         if (err%has_error()) return
         e_active = sum(ci_vectors(:, :, j)*sigma_j)
         gradient_ci_j = 2.0_dp*weights(j)*(sigma_j - e_active*ci_vectors(:, :, j))
         phi = phi + sum(xbar(:, :, j)*gradient_ci_j)
      end do
      call mol%destroy()
   end subroutine ci_response_phi_at

   subroutine solve_z_vector(mol, orbitals, n_inactive, n_active, state, state_index, &
                             kappa_bar, xbar, err)
      !! The genuine Z-vector solution for root `state_index`, exposing
      !! `kappa_bar`/`xbar` directly -- `sa_casscf_gradient_general`'s own RHS
      !! and solve, without the gradient assembly that follows it
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, state_index
      type(sa_hessian_t), intent(in) :: state
      real(dp), allocatable, intent(out) :: kappa_bar(:, :), xbar(:, :, :)
      type(error_t), intent(inout) :: err

      type(mcscf_fock_t) :: fock_i
      real(dp), allocatable :: dm1_i(:, :), dm2_i(:, :, :, :), grad_full(:, :)
      real(dp), allocatable :: rhs(:, :), x(:, :)
      integer :: n_param, l, j, seg0, seg1, na, nb, iterations
      real(dp) :: residual

      if (err%has_error()) return
      call active_space_rdms(state%ci_vectors(:, :, state_index), state%alpha, state%beta, &
                             dm1_i, dm2_i, err)
      if (err%has_error()) return
      call generalized_fock(mol, orbitals, n_inactive, n_active, dm1_i, dm2_i, fock_i, err)
      if (err%has_error()) return
      call orbital_gradient(fock_i, n_inactive, n_active, grad_full)

      n_param = sa_hessian_n_param(state)
      allocate (rhs(n_param, 1))
      rhs = 0.0_dp
      do l = 1, state%n_rot
         rhs(l, 1) = -grad_full(state%rows(l), state%cols(l))
      end do
      call sa_zvector_solve(state, rhs, x, iterations, residual, 1.0e-10_dp, 200, err)
      if (err%has_error()) return

      allocate (kappa_bar(state%n_mo, state%n_mo))
      kappa_bar = 0.0_dp
      do l = 1, state%n_rot
         kappa_bar(state%rows(l), state%cols(l)) = x(l, 1)
         kappa_bar(state%cols(l), state%rows(l)) = -x(l, 1)
      end do

      na = state%alpha%n_strings
      nb = state%beta%n_strings
      allocate (xbar(na, nb, state%n_states))
      do j = 1, state%n_states
         seg0 = state%n_rot + (j - 1)*state%n_det + 1
         seg1 = state%n_rot + j*state%n_det
         xbar(:, :, j) = reshape(x(seg0:seg1, 1), [na, nb])
      end do
   end subroutine solve_z_vector

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
            v(ia, ib) = sin(0.41_dp*real(ia, dp) + 0.29_dp*real(ib, dp) + phase)
         end do
      end do
   end subroutine ci_direction

   ! ============================================================
   ! Gates 1, 2, 5 and the unequal-weights refusal
   ! ============================================================

   subroutine test_short_circuit_bit_identical(error)
      !! `czt_sa_casscf_gradient` at `n_states = 1` must equal
      !! `czt_mcscf_gradient` bit for bit -- it *is* that call
      !!
      !! At one thread: the densities are rebuilt from the CI vector, and a
      !! threaded RDM build's merge order is not fixed from call to call.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), allocatable :: g_sa(:, :), g_base(:, :)
      logical :: ok
      integer :: threads

      threads = omp_get_max_threads()
      call omp_set_num_threads(1)
      call converged_sa(LIH, "sto-3g", 1, 2, [1.0_dp], mol, orbitals, result, err, ok)
      call check(error, ok, "single-state CAS(2,2) should converge")
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if

      call czt_mcscf_gradient(mol, result%orbitals, 1, 2, result%dm1, result%dm2, &
                              g_base, err)
      call check(error,.not. err%has_error(), "czt_mcscf_gradient should not error")
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if

      block
         real(dp), allocatable :: ci3(:, :, :), en1(:)
         allocate (ci3(size(result%ci_vector, 1), size(result%ci_vector, 2), 1))
         ci3(:, :, 1) = result%ci_vector
         en1 = [result%energy]
         call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, ci3, en1, &
                                     [1.0_dp], 1, g_sa, err)
      end block
      call check(error,.not. err%has_error(), "czt_sa_casscf_gradient should not error")
      if (allocated(error)) then
         call omp_set_num_threads(threads)
         return
      end if

      call omp_set_num_threads(threads)
      call check(error, all(g_sa == g_base), &
                 "the n_states=1 short-circuit should be bit-identical")
      if (allocated(error)) return

      ! The short-circuit returns before the general path's range check, so
      ! it carries its own: with one state there is no root 2.
      block
         real(dp), allocatable :: ci3(:, :, :), en1(:), g_bad(:, :)
         type(error_t) :: bad
         allocate (ci3(size(result%ci_vector, 1), size(result%ci_vector, 2), 1))
         ci3(:, :, 1) = result%ci_vector
         en1 = [result%energy]
         call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, ci3, en1, &
                                     [1.0_dp], 2, g_bad, bad)
         call check(error, bad%has_error(), &
                    "root 2 of a one-state average should be refused, not "// &
                    "answered with root 1")
      end block
      call mol%destroy()
   end subroutine test_short_circuit_bit_identical

   subroutine test_general_path_n1(error)
      !! The *general* Z-vector path, forced to run at `n_states = 1` by
      !! calling `sa_casscf_gradient_general` directly (bypassing the literal
      !! short-circuit), should still agree with `czt_mcscf_gradient` -- to
      !! solver tolerance, not bit for bit, since kappa_bar is now ~1e-10
      !! rather than exactly zero
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), allocatable :: g_sa(:, :), g_base(:, :)
      logical :: ok

      call converged_sa(LIH, "sto-3g", 1, 2, [1.0_dp], mol, orbitals, result, err, ok)
      call check(error, ok, "single-state CAS(2,2) should converge")
      if (allocated(error)) return

      call czt_mcscf_gradient(mol, result%orbitals, 1, 2, result%dm1, result%dm2, &
                              g_base, err)
      call check(error,.not. err%has_error(), "czt_mcscf_gradient should not error")
      if (allocated(error)) return

      block
         real(dp), allocatable :: ci3(:, :, :), en1(:)
         allocate (ci3(size(result%ci_vector, 1), size(result%ci_vector, 2), 1))
         ci3(:, :, 1) = result%ci_vector
         en1 = [result%energy]
         call sa_casscf_gradient_general(mol, result%orbitals, 1, 2, 1, 1, ci3, en1, &
                                         [1.0_dp], 1, g_sa, err)
      end block
      call check(error,.not. err%has_error(), "sa_casscf_gradient_general should not error")
      if (allocated(error)) return

      call check(error, maxval(abs(g_sa - g_base)) < 1.0e-8_dp, &
                 "the general Z-vector path at n_states=1 should reproduce "// &
                 "czt_mcscf_gradient to solver tolerance")
      call mol%destroy()
   end subroutine test_general_path_n1

   subroutine test_weighted_sum_identity(error)
      !! `sum_I w_I dE_I/dR = dE_SA/dR` -- see the module docstring for why
      !! this gate is structurally weak
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), allocatable :: g1(:, :), g2(:, :), g_sum(:, :), g_sa(:, :)
      real(dp), allocatable :: dm1_sa(:, :), dm2_sa(:, :, :, :)
      type(link_table_t) :: alpha, beta
      logical :: ok

      call converged_sa(LIH, "sto-3g", 1, 2, SA2_WEIGHTS, mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                  result%energies, SA2_WEIGHTS, 1, g1, err)
      call check(error,.not. err%has_error(), "root 1's gradient should build")
      if (allocated(error)) return
      call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                  result%energies, SA2_WEIGHTS, 2, g2, err)
      call check(error,.not. err%has_error(), "root 2's gradient should build")
      if (allocated(error)) return

      g_sum = SA2_WEIGHTS(1)*g1 + SA2_WEIGHTS(2)*g2

      call build_link_table(2, 1, alpha, err)
      call build_link_table(2, 1, beta, err)
      call check(error,.not. err%has_error(), "the link tables should build")
      if (allocated(error)) return
      call sa_density_matrices(result%ci_vectors, SA2_WEIGHTS, alpha, beta, dm1_sa, &
                               dm2_sa, err)
      call alpha%destroy()
      call beta%destroy()
      call check(error,.not. err%has_error(), "the SA densities should build")
      if (allocated(error)) return
      call czt_mcscf_gradient(mol, result%orbitals, 1, 2, dm1_sa, dm2_sa, g_sa, err)
      call check(error,.not. err%has_error(), "the SA gradient should build")
      if (allocated(error)) return

      call check(error, maxval(abs(g_sum - g_sa)) < 1.0e-8_dp, &
                 "the weighted sum of root gradients should equal the SA gradient")
      call mol%destroy()
   end subroutine test_weighted_sum_identity

   subroutine test_translational_invariance(error)
      !! Necessary, not sufficient: each root's gradient sums to (near) zero
      !! over atoms
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), allocatable :: g1(:, :), g2(:, :)
      logical :: ok
      integer :: comp

      call converged_sa(LIH, "sto-3g", 1, 2, SA2_WEIGHTS, mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                  result%energies, SA2_WEIGHTS, 1, g1, err)
      call check(error,.not. err%has_error(), "root 1's gradient should build")
      if (allocated(error)) return
      call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                  result%energies, SA2_WEIGHTS, 2, g2, err)
      call check(error,.not. err%has_error(), "root 2's gradient should build")
      if (allocated(error)) return

      do comp = 1, 3
         call check(error, abs(sum(g1(comp, :))) < 1.0e-7_dp, &
                    "root 1's gradient should sum to zero over atoms")
         if (allocated(error)) exit
         call check(error, abs(sum(g2(comp, :))) < 1.0e-7_dp, &
                    "root 2's gradient should sum to zero over atoms")
         if (allocated(error)) exit
      end do
      call mol%destroy()
   end subroutine test_translational_invariance

   subroutine test_unequal_weights_refused(error)
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), allocatable :: g(:, :)
      real(dp), parameter :: UNEQUAL(2) = [0.7_dp, 0.3_dp]
      logical :: ok

      call converged_sa(LIH, "sto-3g", 1, 2, UNEQUAL, mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) (0.7, 0.3) should converge")
      if (allocated(error)) return

      call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                  result%energies, UNEQUAL, 1, g, err)
      call check(error, err%has_error(), &
                 "unequal weights should be refused rather than silently accepted")
      call mol%destroy()
   end subroutine test_unequal_weights_refused

   subroutine converged_cas44_sa3(mol, orbitals, result, err, ok)
      !! LiH/6-31G CAS(4,4) SA-3: 18 non-redundant CI directions per state,
      !! the multistate Hessian gate's own system (`test_mqc_sa_hessian_multistate.f90`)
      type(czt_molecule_t), intent(out) :: mol
      real(dp), allocatable, intent(out) :: orbitals(:, :)
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: err
      logical, intent(out) :: ok
      real(dp), parameter :: WEIGHTS3(3) = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]

      type(rhf_result_t) :: scf

      ok = .false.
      call build_czt_molecule(LIH_Z, LIH_SYM, LIH, "6-31g", mol, err)
      if (err%has_error()) return
      call run_czt_rhf(mol, 4, 200, 1.0e-12_dp, 1.0e-10_dp, .false., scf, err)
      if (err%has_error()) return
      if (.not. scf%converged) then
         call err%set(ERROR_VALIDATION, "the LiH/6-31G RHF reference did not converge")
         return
      end if
      call run_czt_casscf(mol, scf%orbitals, 0, 4, 2, 2, result, err, max_iterations=500, &
                          gradient_tol=1.0e-10_dp, n_states=3, weights=WEIGHTS3)
      if (err%has_error()) return
      orbitals = result%orbitals
      ok = result%converged
   end subroutine converged_cas44_sa3

   subroutine test_weighted_sum_cas44_sa3(error)
      !! The weighted-sum identity again, on a system with a genuinely
      !! multi-dimensional CI response (18 directions per state) rather than
      !! the SA-2-CAS(2,2) system's single one -- still structurally weak,
      !! see the module docstring
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: WEIGHTS3(3) = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), allocatable :: g(:, :, :), g_sum(:, :), g_sa(:, :)
      real(dp), allocatable :: dm1_sa(:, :), dm2_sa(:, :, :, :)
      type(link_table_t) :: alpha, beta
      integer :: j
      logical :: ok

      call converged_cas44_sa3(mol, orbitals, result, err, ok)
      call check(error, ok, "CAS(4,4) SA-3 should converge")
      if (allocated(error)) return

      do j = 1, 3
         block
            real(dp), allocatable :: gj(:, :)
            call czt_sa_casscf_gradient(mol, result%orbitals, 0, 4, 2, 2, &
                                        result%ci_vectors, result%energies, WEIGHTS3, j, &
                                        gj, err)
            call check(error,.not. err%has_error(), "each root's gradient should build")
            if (allocated(error)) return
            if (j == 1) then
               allocate (g(3, mol%natm, 3))
            end if
            g(:, :, j) = gj
         end block
      end do

      g_sum = WEIGHTS3(1)*g(:, :, 1) + WEIGHTS3(2)*g(:, :, 2) + WEIGHTS3(3)*g(:, :, 3)

      call build_link_table(4, 2, alpha, err)
      call build_link_table(4, 2, beta, err)
      call check(error,.not. err%has_error(), "the link tables should build")
      if (allocated(error)) return
      call sa_density_matrices(result%ci_vectors, WEIGHTS3, alpha, beta, dm1_sa, dm2_sa, err)
      call alpha%destroy()
      call beta%destroy()
      call check(error,.not. err%has_error(), "the SA densities should build")
      if (allocated(error)) return
      call czt_mcscf_gradient(mol, result%orbitals, 0, 4, dm1_sa, dm2_sa, g_sa, err)
      call check(error,.not. err%has_error(), "the SA gradient should build")
      if (allocated(error)) return

      call check(error, maxval(abs(g_sum - g_sa)) < 1.0e-7_dp, &
                 "the weighted sum of the three root gradients should equal the SA "// &
                 "gradient")
      call mol%destroy()
   end subroutine test_weighted_sum_cas44_sa3

   ! ============================================================
   ! Per-piece gates against the connection-based reference
   ! ============================================================

   subroutine test_base_term_connection(error)
      !! `czt_mcscf_gradient(dm1_I, dm2_I)` against a central difference of
      !! `root_energy_at`, in the connection -- if this fails, the base term
      !! itself is wrong and nothing downstream can be trusted
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 1.0e-4_dp
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: dm1_i(:, :), dm2_i(:, :, :, :), analytic(:, :)
      real(dp) :: moved(3, 2), e_plus, e_minus, fd, worst
      integer :: iatom, comp
      logical :: ok

      call converged_sa(LIH, "sto-3g", 1, 2, SA2_WEIGHTS, mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call build_link_table(2, 1, alpha, err)
      call build_link_table(2, 1, beta, err)
      call active_space_rdms(result%ci_vectors(:, :, 1), alpha, beta, dm1_i, dm2_i, err)
      call alpha%destroy()
      call beta%destroy()
      call check(error,.not. err%has_error(), "root 1's densities should build")
      if (allocated(error)) return

      call czt_mcscf_gradient(mol, result%orbitals, 1, 2, dm1_i, dm2_i, analytic, err)
      call check(error,.not. err%has_error(), "czt_mcscf_gradient should not error")
      if (allocated(error)) return

      worst = 0.0_dp
      do iatom = 1, 2
         do comp = 1, 3
            moved = LIH
            moved(comp, iatom) = LIH(comp, iatom) + STEP
            call root_energy_at(moved, "sto-3g", result%orbitals, 1, 2, dm1_i, dm2_i, &
                                e_plus, err)
            call check(error,.not. err%has_error(), "root_energy_at(+h) should not error")
            if (allocated(error)) exit
            moved(comp, iatom) = LIH(comp, iatom) - STEP
            call root_energy_at(moved, "sto-3g", result%orbitals, 1, 2, dm1_i, dm2_i, &
                                e_minus, err)
            call check(error,.not. err%has_error(), "root_energy_at(-h) should not error")
            if (allocated(error)) exit
            fd = (e_plus - e_minus)/(2.0_dp*STEP)
            worst = max(worst, abs(fd - analytic(comp, iatom)))
         end do
         if (allocated(error)) exit
      end do
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call check(error, worst < 1.0e-7_dp, &
                 "czt_mcscf_gradient should difference root_energy_at in the connection")
      call mol%destroy()
   end subroutine test_base_term_connection

   subroutine test_orbital_response_connection(error)
      !! `orbital_response_gradient` against a central difference of
      !! `orbital_response_phi_at`, both an arbitrary fixed direction and the
      !! genuine Z-vector solution for root 1
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 1.0e-4_dp
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: kappa_bar(:, :), analytic(:, :), xbar_dummy(:, :, :)
      real(dp), allocatable :: kappa_flat(:)
      integer, allocatable :: rows(:), cols(:)
      real(dp) :: moved(3, 2), phi_plus, phi_minus, fd, worst
      integer :: n_rot, l, iatom, comp, kase
      logical :: ok

      call converged_sa(LIH, "sto-3g", 1, 2, SA2_WEIGHTS, mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, SA2_WEIGHTS, state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      n_rot = state%n_rot
      call rotation_parameters(state%n_mo, 1, 2, rows, cols)

      worst = 0.0_dp
      do kase = 1, 2
         if (allocated(kappa_flat)) deallocate (kappa_flat)
         if (kase == 1) then
            kappa_flat = flat_kappa_direction(n_rot, 0.6_dp)
         else
            call solve_z_vector(mol, result%orbitals, 1, 2, state, 1, kappa_bar, &
                                xbar_dummy, err)
            call check(error,.not. err%has_error(), "solve_z_vector should not error")
            if (allocated(error)) exit
            allocate (kappa_flat(n_rot))
            do l = 1, n_rot
               kappa_flat(l) = kappa_bar(rows(l), cols(l))
            end do
            deallocate (kappa_bar)
         end if

         if (allocated(kappa_bar)) deallocate (kappa_bar)
         allocate (kappa_bar(state%n_mo, state%n_mo))
         kappa_bar = 0.0_dp
         do l = 1, n_rot
            kappa_bar(rows(l), cols(l)) = kappa_flat(l)
            kappa_bar(cols(l), rows(l)) = -kappa_flat(l)
         end do

         call orbital_response_gradient(mol, result%orbitals, 1, 2, state, kappa_bar, &
                                        analytic, err)
         call check(error,.not. err%has_error(), "orbital_response_gradient should not error")
         if (allocated(error)) exit

         do iatom = 1, 2
            do comp = 1, 3
               moved = LIH
               moved(comp, iatom) = LIH(comp, iatom) + STEP
               call orbital_response_phi_at(moved, "sto-3g", result%orbitals, 1, 2, &
                                            state%dm1_sa, state%dm2_sa, kappa_flat, rows, &
                                            cols, phi_plus, err)
               call check(error,.not. err%has_error(), "phi_at(+h) should not error")
               if (allocated(error)) exit
               moved(comp, iatom) = LIH(comp, iatom) - STEP
               call orbital_response_phi_at(moved, "sto-3g", result%orbitals, 1, 2, &
                                            state%dm1_sa, state%dm2_sa, kappa_flat, rows, &
                                            cols, phi_minus, err)
               call check(error,.not. err%has_error(), "phi_at(-h) should not error")
               if (allocated(error)) exit
               fd = (phi_plus - phi_minus)/(2.0_dp*STEP)
               worst = max(worst, abs(fd - analytic(comp, iatom)))
            end do
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
         deallocate (kappa_bar)
      end do
      if (allocated(error)) then
         call destroy_sa_hessian(state)
         call mol%destroy()
         return
      end if

      call check(error, worst < 1.0e-6_dp, &
                 "orbital_response_gradient should difference orbital_response_phi_at "// &
                 "in the connection")
      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_orbital_response_connection

   subroutine test_ci_response_connection(error)
      !! `ci_response_gradient` against a central difference of
      !! `ci_response_phi_at`, both an arbitrary fixed (projected, symmetric)
      !! direction and the genuine Z-vector solution for root 1
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 1.0e-4_dp
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      type(sa_hessian_t) :: state
      real(dp), allocatable :: xbar(:, :, :), analytic(:, :), kappa_bar_dummy(:, :)
      real(dp) :: moved(3, 2), phi_plus, phi_minus, fd, worst
      integer :: na, nb, j, iatom, comp, kase
      logical :: ok

      call converged_sa(LIH, "sto-3g", 1, 2, SA2_WEIGHTS, mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call build_sa_hessian(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                            result%energies, SA2_WEIGHTS, state, err)
      call check(error,.not. err%has_error(), "build_sa_hessian should not error")
      if (allocated(error)) return

      na = state%alpha%n_strings
      nb = state%beta%n_strings

      worst = 0.0_dp
      do kase = 1, 2
         if (allocated(xbar)) deallocate (xbar)
         if (kase == 1) then
            allocate (xbar(na, nb, 2))
            do j = 1, 2
               call ci_direction(na, nb, 0.8_dp + 0.2_dp*real(j, dp), xbar(:, :, j))
               call project_ci_block(state, xbar(:, :, j))
               xbar(:, :, j) = 0.5_dp*(xbar(:, :, j) + transpose(xbar(:, :, j)))
            end do
         else
            call solve_z_vector(mol, result%orbitals, 1, 2, state, 1, kappa_bar_dummy, &
                                xbar, err)
            call check(error,.not. err%has_error(), "solve_z_vector should not error")
            if (allocated(error)) exit
            deallocate (kappa_bar_dummy)
         end if

         ! The Z-vector solution carries components along the averaged states
         ! (null directions of the projected Hessian). Their Lagrange
         ! multipliers are zero for non-degenerate states, so the reference is
         ! taken with the projected vector, as `ci_response_gradient` uses it.
         ! Left unprojected, the reference picks up d/dR <c_K|H|c_J>, which is
         ! not part of the gradient.
         if (kase == 2) then
            do j = 1, 2
               call project_ci_block(state, xbar(:, :, j))
            end do
         end if

         call ci_response_gradient(mol, result%orbitals, 1, 2, state, xbar, SA2_WEIGHTS, &
                                   analytic, err)
         call check(error,.not. err%has_error(), "ci_response_gradient should not error")
         if (allocated(error)) exit

         do iatom = 1, 2
            do comp = 1, 3
               moved = LIH
               moved(comp, iatom) = LIH(comp, iatom) + STEP
               call ci_response_phi_at(moved, "sto-3g", result%orbitals, 1, 2, 1, 1, &
                                       result%ci_vectors, SA2_WEIGHTS, xbar, state%alpha, &
                                       state%beta, phi_plus, err)
               call check(error,.not. err%has_error(), "phi_at(+h) should not error")
               if (allocated(error)) exit
               moved(comp, iatom) = LIH(comp, iatom) - STEP
               call ci_response_phi_at(moved, "sto-3g", result%orbitals, 1, 2, 1, 1, &
                                       result%ci_vectors, SA2_WEIGHTS, xbar, state%alpha, &
                                       state%beta, phi_minus, err)
               call check(error,.not. err%has_error(), "phi_at(-h) should not error")
               if (allocated(error)) exit
               fd = (phi_plus - phi_minus)/(2.0_dp*STEP)
               worst = max(worst, abs(fd - analytic(comp, iatom)))
            end do
            if (allocated(error)) exit
         end do
         if (allocated(error)) exit
      end do
      if (allocated(error)) then
         call destroy_sa_hessian(state)
         call mol%destroy()
         return
      end if

      call check(error, worst < 1.0e-6_dp, &
                 "ci_response_gradient should difference ci_response_phi_at in the "// &
                 "connection")
      call destroy_sa_hessian(state)
      call mol%destroy()
   end subroutine test_ci_response_connection

   ! ============================================================
   ! Gate 3: finite differences of this code's own re-converged energies
   ! ============================================================

   subroutine energy_at_lih(coordinates, basis, n_inactive, n_active, weights, energies, err)
      !! Every root's energy at a displaced LiH geometry, a fresh RHF/CASSCF
      real(dp), intent(in) :: coordinates(3, 2)
      character(len=*), intent(in) :: basis
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: weights(:)
      real(dp), allocatable, intent(out) :: energies(:)
      type(error_t), intent(inout) :: err

      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      logical :: ok

      call converged_sa(coordinates, basis, n_inactive, n_active, weights, mol, orbitals, &
                        result, err, ok)
      if (err%has_error()) return
      if (.not. ok) then
         call err%set(ERROR_VALIDATION, "the displaced-geometry SA-CASSCF did not converge")
         return
      end if
      energies = result%energies
      call mol%destroy()
   end subroutine energy_at_lih

   subroutine test_finite_difference_lih_sa2(error)
      !! Central differences of `run_czt_casscf`'s own root energies against
      !! the analytic per-root gradient, every coordinate, LiH/STO-3G SA-2
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 1.0e-3_dp
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), allocatable :: g1(:, :), g2(:, :)
      real(dp), allocatable :: e_plus(:), e_minus(:)
      real(dp) :: moved(3, 2), fd1, fd2, worst1, worst2
      integer :: iatom, comp
      logical :: ok

      call converged_sa(LIH, "sto-3g", 1, 2, SA2_WEIGHTS, mol, orbitals, result, err, ok)
      call check(error, ok, "SA-2-CAS(2,2) should converge")
      if (allocated(error)) return

      call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                  result%energies, SA2_WEIGHTS, 1, g1, err)
      call check(error,.not. err%has_error(), "root 1's gradient should build")
      if (allocated(error)) return
      call czt_sa_casscf_gradient(mol, result%orbitals, 1, 2, 1, 1, result%ci_vectors, &
                                  result%energies, SA2_WEIGHTS, 2, g2, err)
      call check(error,.not. err%has_error(), "root 2's gradient should build")
      if (allocated(error)) return

      worst1 = 0.0_dp
      worst2 = 0.0_dp
      outer: do iatom = 1, 2
         do comp = 1, 3
            moved = LIH
            moved(comp, iatom) = LIH(comp, iatom) + STEP
            call energy_at_lih(moved, "sto-3g", 1, 2, SA2_WEIGHTS, e_plus, err)
            call check(error,.not. err%has_error(), "the +h displacement should converge")
            if (allocated(error)) exit outer
            moved(comp, iatom) = LIH(comp, iatom) - STEP
            call energy_at_lih(moved, "sto-3g", 1, 2, SA2_WEIGHTS, e_minus, err)
            call check(error,.not. err%has_error(), "the -h displacement should converge")
            if (allocated(error)) exit outer

            fd1 = (e_plus(1) - e_minus(1))/(2.0_dp*STEP)
            fd2 = (e_plus(2) - e_minus(2))/(2.0_dp*STEP)
            worst1 = max(worst1, abs(fd1 - g1(comp, iatom)))
            worst2 = max(worst2, abs(fd2 - g2(comp, iatom)))
         end do
      end do outer
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call check(error, worst1 < 1.0e-6_dp, &
                 "root 1's analytic gradient should difference its own energy")
      if (allocated(error)) then
         call mol%destroy()
         return
      end if
      call check(error, worst2 < 1.0e-6_dp, &
                 "root 2's analytic gradient should difference its own energy")
      call mol%destroy()
   end subroutine test_finite_difference_lih_sa2

   subroutine test_finite_difference_cas44_sa3(error)
      !! Finite differences on the bond-length coordinate of the CAS(4,4)
      !! SA-3 system, every root, so the CI-response and orbital-CI blocks
      !! (18 directions per state) are exercised, not only the single
      !! direction SA-2-CAS(2,2) has
      type(error_type), allocatable, intent(out) :: error
      real(dp), parameter :: STEP = 1.0e-3_dp
      real(dp), parameter :: WEIGHTS3(3) = [1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
      type(error_t) :: err
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)
      type(casscf_result_t) :: result
      real(dp), allocatable :: g(:, :, :), g_root(:, :)
      real(dp), allocatable :: e_plus(:), e_minus(:)
      real(dp) :: moved(3, 2), fd, worst
      integer :: iatom, comp, root
      logical :: ok

      call converged_cas44_sa3(mol, orbitals, result, err, ok)
      call check(error, ok, "CAS(4,4) SA-3 should converge")
      if (allocated(error)) return

      allocate (g(3, 2, 3))
      do root = 1, 3
         call czt_sa_casscf_gradient(mol, result%orbitals, 0, 4, 2, 2, result%ci_vectors, &
                                     result%energies, WEIGHTS3, root, g_root, err)
         call check(error,.not. err%has_error(), "every root's gradient should build")
         if (allocated(error)) return
         g(:, :, root) = g_root
      end do

      worst = 0.0_dp
      do iatom = 1, 2
         comp = 3
         moved = LIH
         moved(comp, iatom) = LIH(comp, iatom) + STEP
         call energy_at_lih(moved, "6-31g", 0, 4, WEIGHTS3, e_plus, err)
         call check(error,.not. err%has_error(), "the +h displacement should converge")
         if (allocated(error)) exit
         moved(comp, iatom) = LIH(comp, iatom) - STEP
         call energy_at_lih(moved, "6-31g", 0, 4, WEIGHTS3, e_minus, err)
         call check(error,.not. err%has_error(), "the -h displacement should converge")
         if (allocated(error)) exit

         do root = 1, 3
            fd = (e_plus(root) - e_minus(root))/(2.0_dp*STEP)
            worst = max(worst, abs(fd - g(comp, iatom, root)))
         end do
      end do
      if (allocated(error)) then
         call mol%destroy()
         return
      end if

      call check(error, worst < 1.0e-6_dp, &
                 "every root's analytic gradient should difference its own energy "// &
                 "on the CAS(4,4) SA-3 system")
      call mol%destroy()
   end subroutine test_finite_difference_cas44_sa3

end module test_mqc_sa_gradient

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_sa_gradient, only: collect_mqc_sa_gradient_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_sa_gradient", collect_mqc_sa_gradient_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
