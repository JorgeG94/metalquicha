!! The analytic nuclear gradient of one root of a state-averaged CASSCF
module mqc_czt_sa_gradient
   !!
   !! For root `I` of an SA-CASSCF converged by `run_czt_casscf`, the gradient
   !! is the explicit R-derivative of the Lagrangian
   !!
   !!     L_I = E_I + kappa_bar . dE_SA/dkappa + sum_J xbar_J . g_CI,J
   !!
   !! with `(kappa_bar, xbar)` from the Z-vector equation `H_SA z = -dE_I/dz`
   !! (`sa_zvector_solve`, on `mqc_czt_sa_hessian`'s operator). It is the sum
   !! of three pieces:
   !!
   !! 1. `czt_mcscf_gradient` on root `I`'s own densities at the SA orbitals;
   !! 2. `orbital_response_gradient`, from `kappa_bar`;
   !! 3. `ci_response_gradient`, from `xbar` projected off every averaged state.
   !!
   !! Each piece has its own overlap term. The orbitals follow the symmetric
   !! connection `C(R) = C0 (C0^T S(R) C0)^(-1/2)`, the same one behind
   !! `czt_mcscf_gradient`'s Pulay term. The derivation is in
   !! `mqc_docs/source/developer_sa_casscf.rst`, "Relaxed density assembly".
   !! Equal weights only: unequal ones are refused by name, because the
   !! projection inside `H_SA` is exact only for equal weights.
   !!
   !! `czt_sa_casscf_gradients` computes several roots together. It builds the
   !! SA Hessian state once and solves every root's Z-vector equation in one
   !! block PCG, in which a converged column stops costing Hessian-vector
   !! products. It then shares the two expensive derivative-integral sweeps
   !! across every root and every response piece, rather than repeating them
   !! once per root: one `two_electron_deriv_many` call over the union of
   !! every root's separable densities (`response_separable_assemble`,
   !! `base_gradient_assemble` do the per-root/per-piece dot products from its
   !! potentials), and one stacked sweep over the active two-body Gamma
   !! tensors (`active_two_electron_gradient_stacked`, on
   !! `active_two_electron_gradient_many`). The one-electron core-Hamiltonian
   !! derivative (`build_core_hamiltonian_derivative`) is shared the same way.
   !! `sa_casscf_gradient_general` (`czt_sa_casscf_gradient`'s single-root
   !! path) does not share any of this: its response pieces
   !! (`orbital_response_gradient`, `ci_response_gradient`) each run their own
   !! `two_electron_deriv_many` sweep, which the two `response_separable_*`
   !! routines below only refactor apart, not fuse -- the fusion is
   !! `sa_casscf_gradients_general`'s alone.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks, shell_dim
   use mqc_czt_gradient, only: one_electron_deriv, iprinv_deriv_at, two_electron_deriv_many, &
                               nuclear_repulsion_gradient, DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_czt_mcscf, only: mcscf_fock_t, generalized_fock, orbital_gradient, one_index_fock
   use mqc_czt_mcscf_gradient, only: czt_mcscf_gradient, active_two_electron_gradient, &
                                     active_two_electron_gradient_response, gamma_block, &
                                     active_two_electron_gradient_many, &
                                     cumulant_two_particle_density
   use mqc_czt_sa_hessian, only: sa_hessian_t, build_sa_hessian, destroy_sa_hessian, &
                                 sa_hessian_n_param, sa_hessian_apply, &
                                 sa_hessian_precondition, cheap_generalized_fock, &
                                 project_ci_block
   use mqc_determinants, only: build_link_table, link_table_t
   use mqc_rdm, only: active_space_rdms, transition_rdms
   implicit none
   private

   public :: czt_sa_casscf_gradient
   public :: sa_casscf_gradient_general
   !! Exposed for the tests: the Z-vector path `czt_sa_casscf_gradient`
      !! takes for `n_states > 1`, callable directly at `n_states = 1` (with
      !! `state_index = 1`) to check that the *general* machinery reproduces
      !! `czt_mcscf_gradient` too, not only the literal short-circuit.
   public :: sa_zvector_solve   !! Exposed for the tests
   public :: response_separable_gradient   !! Exposed for the tests
   public :: orbital_response_gradient   !! Exposed for the tests
   public :: ci_response_gradient   !! Exposed for the tests

   public :: czt_sa_casscf_gradients
      !! Every root in `roots`, computed together (one SA Hessian state, one
      !! block Z-vector solve). Per root, agrees with `czt_sa_casscf_gradient`
      !! to the solver tolerance.
   public :: sa_casscf_gradients_general   !! Exposed for the tests
   public :: sa_block_zvector_solve   !! Exposed for the tests
   public :: sa_gradients_on_state
      !! Roots plus caller-built Lagrangian columns on one SA Hessian state;
      !! `mqc_czt_sa_nac` fuses its pairs in through this.
   public :: UNEQUAL_WEIGHT_TOL
      !! Exposed so a caller (`mqc_czt_bridge`) can refuse unequal weights
      !! early, before running the orbital optimisation at all, with the same
      !! threshold this module refuses them with internally.

   real(dp), parameter :: DEFAULT_CG_TOL = 1.0e-10_dp
   integer, parameter :: DEFAULT_CG_MAX_ITER = 200
   real(dp), parameter :: UNEQUAL_WEIGHT_TOL = 1.0e-8_dp

contains

   subroutine czt_sa_casscf_gradient(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                     ci_vectors, energies, weights, state_index, &
                                     gradient, error, cg_tol, cg_max_iter, cg_iterations, &
                                     cg_residual)
      !! The analytic gradient of root `state_index` of a converged SA-CASSCF
      !!
      !! `n_states = 1` (`size(weights) == 1`) is a literal short-circuit to
      !! `czt_mcscf_gradient`, bit for bit -- gate 1 of the phase-4 plan.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)         !! (n_ao, n_mo), the SA orbitals
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)    !! (na, nb, >= n_states)
      real(dp), intent(in) :: energies(:)            !! (>= n_states), total
      real(dp), intent(in) :: weights(:)              !! (n_states)
      integer, intent(in) :: state_index              !! Which root, 1-based
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, n_atoms)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: cg_tol
         !! Relative residual the Z-vector solve stops at. Default `1e-10`.
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations
      real(dp), intent(out), optional :: cg_residual
         !! The relative residual actually reached

      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: dm1_i(:, :), dm2_i(:, :, :, :)

      if (error%has_error()) return
      if (present(cg_iterations)) cg_iterations = 0
      if (present(cg_residual)) cg_residual = 0.0_dp

      if (size(weights) == 1) then
         call build_link_table(n_active, n_alpha, alpha, error)
         call build_link_table(n_active, n_beta, beta, error)
         if (error%has_error()) return
         call active_space_rdms(ci_vectors(:, :, 1), alpha, beta, dm1_i, dm2_i, error)
         call alpha%destroy()
         call beta%destroy()
         if (error%has_error()) return
         call czt_mcscf_gradient(mol, orbitals, n_inactive, n_active, dm1_i, dm2_i, &
                                 gradient, error)
         return
      end if

      call sa_casscf_gradient_general(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                      ci_vectors, energies, weights, state_index, gradient, &
                                      error, cg_tol, cg_max_iter, cg_iterations, cg_residual)
   end subroutine czt_sa_casscf_gradient

   subroutine sa_casscf_gradient_general(mol, orbitals, n_inactive, n_active, n_alpha, &
                                         n_beta, ci_vectors, energies, weights, state_index, &
                                         gradient, error, cg_tol, cg_max_iter, cg_iterations, &
                                         cg_residual)
      !! The Z-vector machinery itself, with no `n_states = 1` short-circuit --
      !! see `czt_sa_casscf_gradient`, which is this with that short-circuit in
      !! front of it. Callable directly at `n_states = 1` to check that the
      !! *general* path also reproduces `czt_mcscf_gradient` (to solver
      !! tolerance, not bit for bit -- gate 1's "separately" clause).
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)
      real(dp), intent(in) :: energies(:)
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: state_index
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: cg_tol
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations
      real(dp), intent(out), optional :: cg_residual

      type(sa_hessian_t) :: state
      real(dp), allocatable :: dm1_i(:, :), dm2_i(:, :, :, :)
      type(mcscf_fock_t) :: fock_i
      real(dp), allocatable :: grad_full(:, :)
      real(dp), allocatable :: rhs(:, :), x(:, :)
      real(dp), allocatable :: kappa_bar(:, :), xbar(:, :, :)
      real(dp), allocatable :: base_gradient(:, :), orb_gradient(:, :), ci_gradient(:, :)
      real(dp) :: use_tol, residual
      integer :: n_states, use_max_iter, iterations
      integer :: n_mo, na, nb, n_param, l, j, seg0, seg1

      if (error%has_error()) return
      n_states = size(weights)
      if (present(cg_iterations)) cg_iterations = 0
      if (present(cg_residual)) cg_residual = 0.0_dp

      if (state_index < 1 .or. state_index > n_states) then
         call error%set(ERROR_VALIDATION, "sa_casscf_gradient: state "// &
                        to_char(state_index)//" is not one of the "// &
                        to_char(n_states)//" averaged states.")
         return
      end if
      if (maxval(weights) - minval(weights) > UNEQUAL_WEIGHT_TOL) then
         call error%set(ERROR_VALIDATION, "sa_casscf_gradient: unequal SA weights -- "// &
                        "the redundancy projection in the SA Hessian is exact only for "// &
                        "equal weights, so a single root's gradient is refused rather "// &
                        "than built from the wrong Lagrangian.")
         return
      end if

      use_tol = DEFAULT_CG_TOL
      if (present(cg_tol)) use_tol = cg_tol
      use_max_iter = DEFAULT_CG_MAX_ITER
      if (present(cg_max_iter)) use_max_iter = cg_max_iter

      n_mo = size(orbitals, 2)

      call build_sa_hessian(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                            ci_vectors, energies, weights, state, error)
      if (error%has_error()) return

      ! ---- the right-hand side: root I's own orbital gradient at the SA
      ! orbitals, and an identically-zero CI part ----------------------------
      call active_space_rdms(ci_vectors(:, :, state_index), state%alpha, state%beta, &
                             dm1_i, dm2_i, error)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if
      call generalized_fock(mol, orbitals, n_inactive, n_active, dm1_i, dm2_i, fock_i, error)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if
      call orbital_gradient(fock_i, n_inactive, n_active, grad_full)

      n_param = sa_hessian_n_param(state)
      allocate (rhs(n_param, 1))
      rhs = 0.0_dp
      do l = 1, state%n_rot
         rhs(l, 1) = -grad_full(state%rows(l), state%cols(l))
      end do

      call sa_zvector_solve(state, rhs, x, iterations, residual, use_tol, use_max_iter, error)
      if (present(cg_iterations)) cg_iterations = iterations
      if (present(cg_residual)) cg_residual = residual
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if

      allocate (kappa_bar(n_mo, n_mo))
      kappa_bar = 0.0_dp
      do l = 1, state%n_rot
         kappa_bar(state%rows(l), state%cols(l)) = x(l, 1)
         kappa_bar(state%cols(l), state%rows(l)) = -x(l, 1)
      end do

      na = state%alpha%n_strings
      nb = state%beta%n_strings
      allocate (xbar(na, nb, n_states))
      do j = 1, n_states
         seg0 = state%n_rot + (j - 1)*state%n_det + 1
         seg1 = state%n_rot + j*state%n_det
         xbar(:, :, j) = reshape(x(seg0:seg1, 1), [na, nb])
      end do

      call czt_mcscf_gradient(mol, orbitals, n_inactive, n_active, dm1_i, dm2_i, &
                              base_gradient, error)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if

      call orbital_response_gradient(mol, orbitals, n_inactive, n_active, state, kappa_bar, &
                                     orb_gradient, error)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if

      call ci_response_gradient(mol, orbitals, n_inactive, n_active, state, xbar, weights, &
                                ci_gradient, error)
      if (error%has_error()) then
         call destroy_sa_hessian(state)
         return
      end if

      gradient = base_gradient + orb_gradient + ci_gradient

      call destroy_sa_hessian(state)
   end subroutine sa_casscf_gradient_general

   subroutine sa_zvector_solve(state, rhs, x, iterations, residual, tol, max_iter, error)
      !! Preconditioned CG on `H_SA x = rhs` for a single right-hand side:
      !! `sa_block_zvector_solve` with one column
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: rhs(:, :)      !! (n_param, 1)
      real(dp), allocatable, intent(out) :: x(:, :)   !! (n_param, 1)
      integer, intent(out) :: iterations
      real(dp), intent(out) :: residual   !! The relative residual reached
      real(dp), intent(in) :: tol
      integer, intent(in) :: max_iter
      type(error_t), intent(inout) :: error

      integer :: its(1)
      real(dp) :: res(1)

      iterations = 0
      residual = 0.0_dp
      if (error%has_error()) return
      if (size(rhs, 2) /= 1) then
         call error%set(ERROR_VALIDATION, "sa_zvector_solve takes one right-hand side; "// &
                        "use sa_block_zvector_solve for "//to_char(size(rhs, 2)))
         return
      end if
      call sa_block_zvector_solve(state, rhs, x, its, res, tol, max_iter, error)
      iterations = its(1)
      residual = res(1)
   end subroutine sa_zvector_solve

   subroutine sa_block_zvector_solve(state, rhs, x, iterations, residual, tol, max_iter, error)
      !! `H_SA x_J = rhs_J` for every column `J` of `rhs`, solved together but
      !! **as `n_vec` independent CG recursions**, not as one CG on the
      !! concatenated vector (which is what calling `sa_hessian_apply` with a
      !! multi-column block and reducing with a single scalar `sum(r*z)` would
      !! do, and what `sa_zvector_solve` itself does when handed `n_vec > 1` --
      !! it is exposed with `n_vec = 1` in every existing call, so that path is
      !! untouched). The SA Hessian does not couple different roots'
      !! right-hand sides, so each column may take its own step size and its
      !! own conjugate-direction history, and this reduces to `sa_zvector_solve`
      !! called once per column to the solver tolerance, at a fraction of the
      !! cost: **a column stops costing Hessian-vector products once its own
      !! relative residual is below `tol`**, tracked per column and gathered
      !! into a shrinking active set before every `sa_hessian_apply` call.
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: rhs(:, :)         !! (n_param, n_vec)
      real(dp), allocatable, intent(out) :: x(:, :)
      integer, intent(out) :: iterations(:)     !! (n_vec), the iteration each column stopped at
      real(dp), intent(out) :: residual(:)      !! (n_vec), the relative residual each reached
      real(dp), intent(in) :: tol
      integer, intent(in) :: max_iter
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: r(:, :), z(:, :), p(:, :), ap(:, :)
      real(dp), allocatable :: p_live(:, :), ap_live(:, :), r_live(:, :), z_live(:, :)
      real(dp), allocatable :: rz(:), rz_new(:), target_norm(:)
      logical, allocatable :: active(:)
      integer, allocatable :: idx(:)
      real(dp) :: pap, step, r_norm
      integer :: n_param, n_vec, iter, iv, j, n_live

      if (error%has_error()) return
      n_param = size(rhs, 1)
      n_vec = size(rhs, 2)
      allocate (x(n_param, n_vec), r(n_param, n_vec), z(n_param, n_vec))
      allocate (p(n_param, n_vec), ap(n_param, n_vec))
      allocate (rz(n_vec), rz_new(n_vec), target_norm(n_vec), active(n_vec))

      x = 0.0_dp
      r = rhs
      iterations = 0
      residual = 0.0_dp
      do iv = 1, n_vec
         target_norm(iv) = sqrt(sum(rhs(:, iv)**2))
      end do
      active = target_norm > 0.0_dp
      if (.not. any(active)) return

      call sa_hessian_precondition(state, r, z)
      p = z
      do iv = 1, n_vec
         rz(iv) = sum(r(:, iv)*z(:, iv))
      end do

      do iter = 1, max_iter
         idx = pack([(iv, iv=1, n_vec)], active)
         n_live = size(idx)
         if (n_live == 0) exit
         do j = 1, n_live
            iterations(idx(j)) = iter
         end do

         ! A vector-subscripted section (`p(:, idx)`) cannot itself be an
         ! INTENT(OUT)/INTENT(INOUT) actual argument, so the active columns
         ! are gathered into a contiguous scratch array and scattered back --
         ! copies `sa_hessian_apply`'s own per-column loop does not need,
         ! since it is the caller here that is discontiguous, not the callee.
         allocate (p_live(n_param, n_live), ap_live(n_param, n_live))
         do j = 1, n_live
            p_live(:, j) = p(:, idx(j))
         end do
         call sa_hessian_apply(state, p_live, ap_live, error)
         do j = 1, n_live
            ap(:, idx(j)) = ap_live(:, j)
         end do
         deallocate (p_live, ap_live)
         if (error%has_error()) return

         do j = 1, n_live
            iv = idx(j)
            pap = sum(p(:, iv)*ap(:, iv))
            if (pap <= 0.0_dp) then
               call error%set(ERROR_VALIDATION, "sa_casscf_gradients: the Z-vector CG "// &
                              "search direction has non-positive curvature, so the SA "// &
                              "Hessian is not positive definite on the non-redundant "// &
                              "space -- the reference is not a genuine SA-CASSCF minimum.")
               return
            end if
            step = rz(iv)/pap
            x(:, iv) = x(:, iv) + step*p(:, iv)
            r(:, iv) = r(:, iv) - step*ap(:, iv)
            r_norm = sqrt(sum(r(:, iv)**2))
            residual(iv) = r_norm/target_norm(iv)
            if (residual(iv) <= tol) active(iv) = .false.
         end do

         if (.not. any(active)) exit
         idx = pack([(iv, iv=1, n_vec)], active)
         n_live = size(idx)
         allocate (r_live(n_param, n_live), z_live(n_param, n_live))
         do j = 1, n_live
            r_live(:, j) = r(:, idx(j))
         end do
         call sa_hessian_precondition(state, r_live, z_live)
         do j = 1, n_live
            z(:, idx(j)) = z_live(:, j)
         end do
         deallocate (r_live, z_live)
         do j = 1, n_live
            iv = idx(j)
            rz_new(iv) = sum(r(:, iv)*z(:, iv))
            p(:, iv) = z(:, iv) + (rz_new(iv)/rz(iv))*p(:, iv)
            rz(iv) = rz_new(iv)
         end do
      end do

      deallocate (r, z, p, ap, rz, rz_new, target_norm, active)
   end subroutine sa_block_zvector_solve

   subroutine czt_sa_casscf_gradients(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                      ci_vectors, energies, weights, roots, gradients, error, &
                                      cg_tol, cg_max_iter, cg_iterations, cg_residual)
      !! The analytic gradients of every root named in `roots`, computed
      !! together: `build_sa_hessian` once (rather than once per root, as
      !! `n_roots` separate `czt_sa_casscf_gradient` calls would), the
      !! Z-vector solved for every root's right-hand side in one adaptive
      !! block PCG (`sa_block_zvector_solve`), and each root's relaxed
      !! gradient assembled by `sa_casscf_gradients_general`, which shares
      !! every derivative-integral sweep across every root -- see that
      !! routine's own docstring for what is shared and how.
      !!
      !! `n_states = 1` (`size(weights) == 1`) is the same bit-identical
      !! short-circuit to `czt_mcscf_gradient` `czt_sa_casscf_gradient` takes,
      !! applied to whichever root(s) of `roots` are asked for (there is only
      !! one root to ask for).
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)         !! (n_ao, n_mo), the SA orbitals
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)    !! (na, nb, >= n_states)
      real(dp), intent(in) :: energies(:)            !! (>= n_states), total
      real(dp), intent(in) :: weights(:)              !! (n_states)
      integer, intent(in) :: roots(:)                  !! 1-based, the roots asked for
      real(dp), allocatable, intent(out) :: gradients(:, :, :)  !! (3, n_atoms, size(roots))
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: cg_tol
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations(:)   !! (size(roots))
      real(dp), intent(out), optional :: cg_residual(:)     !! (size(roots))

      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: dm1_i(:, :), dm2_i(:, :, :, :)
      real(dp), allocatable :: gradient(:, :)
      integer :: ir

      if (error%has_error()) return

      if (size(weights) == 1) then
         call build_link_table(n_active, n_alpha, alpha, error)
         call build_link_table(n_active, n_beta, beta, error)
         if (error%has_error()) return
         call active_space_rdms(ci_vectors(:, :, 1), alpha, beta, dm1_i, dm2_i, error)
         call alpha%destroy()
         call beta%destroy()
         if (error%has_error()) return
         allocate (gradients(3, mol%natm, size(roots)))
         do ir = 1, size(roots)
            call czt_mcscf_gradient(mol, orbitals, n_inactive, n_active, dm1_i, dm2_i, &
                                    gradient, error)
            if (error%has_error()) return
            gradients(:, :, ir) = gradient
         end do
         return
      end if

      call sa_casscf_gradients_general(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                       ci_vectors, energies, weights, roots, gradients, error, &
                                       cg_tol, cg_max_iter, cg_iterations, cg_residual)
   end subroutine czt_sa_casscf_gradients

   subroutine sa_casscf_gradients_general(mol, orbitals, n_inactive, n_active, n_alpha, &
                                          n_beta, ci_vectors, energies, weights, roots, &
                                          gradients, error, cg_tol, cg_max_iter, &
                                          cg_iterations, cg_residual)
      !! The fused Z-vector machinery `czt_sa_casscf_gradients` takes for
      !! `n_states > 1`: `build_sa_hessian`, then `sa_gradients_on_state`
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)
      real(dp), intent(in) :: energies(:)
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: roots(:)
      real(dp), allocatable, intent(out) :: gradients(:, :, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: cg_tol
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations(:)
      real(dp), intent(out), optional :: cg_residual(:)

      type(sa_hessian_t) :: state

      if (error%has_error()) return
      call check_sa_request(size(weights), weights, roots, error)
      if (error%has_error()) return
      call build_sa_hessian(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                            ci_vectors, energies, weights, state, error)
      if (error%has_error()) return
      call sa_gradients_on_state(state, orbitals, n_inactive, n_active, weights, roots, &
                                 gradients, error, cg_tol, cg_max_iter, cg_iterations, &
                                 cg_residual)
      call destroy_sa_hessian(state)
   end subroutine sa_casscf_gradients_general

   subroutine check_sa_request(n_states, weights, roots, error)
      !! Refuse a root outside the averaged states, and unequal weights
      integer, intent(in) :: n_states
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: roots(:)
      type(error_t), intent(inout) :: error

      integer :: ir

      do ir = 1, size(roots)
         if (roots(ir) < 1 .or. roots(ir) > n_states) then
            call error%set(ERROR_VALIDATION, "sa_casscf_gradients: root "// &
                           to_char(roots(ir))//" is not one of the "// &
                           to_char(n_states)//" averaged states.")
            return
         end if
      end do
      if (maxval(weights) - minval(weights) > UNEQUAL_WEIGHT_TOL) then
         call error%set(ERROR_VALIDATION, "sa_casscf_gradients: unequal SA weights -- "// &
                        "the redundancy projection in the SA Hessian is exact only for "// &
                        "equal weights, so a single root's gradient is refused rather "// &
                        "than built from the wrong Lagrangian.")
      end if
   end subroutine check_sa_request

   subroutine sa_gradients_on_state(state, orbitals, n_inactive, n_active, weights, roots, &
                                    gradients, error, cg_tol, cg_max_iter, cg_iterations, &
                                    cg_residual, extra_rhs, extra_gamma, extra_d_active, &
                                    extra_weighted, extra_out)
      !! Every root in `roots`, and any extra Lagrangian columns, relaxed
      !! together on an already-built SA Hessian `state`
      !!
      !! A root column is `L_I`; its right-hand side is built here. An extra
      !! column is a Lagrangian whose right-hand side (`extra_rhs`) and base
      !! term the caller supplies: an active two-body density on `c_active`
      !! (`extra_gamma`), and an active one-body AO density with its
      !! energy-weighted matrix (`extra_d_active`, `extra_weighted`) in the
      !! same form as the CI-response piece. It has no nuclear-repulsion or
      !! core-density base; a NAC pair's `<I|dH/dR|J>` is exactly this.
      !! `extra_out` is each extra column's assembled derivative.
      !!
      !! **What is shared.** One block PCG over every column, in which a
      !! converged column stops costing Hessian-vector products. Every
      !! column's cheap (no derivative-integral) densities and energy-weighted
      !! matrices are built first (`orbital_response_pieces`,
      !! `ci_response_pieces`, `cumulant_two_particle_density`,
      !! `build_active_density`, `build_weighted_from_fock`). Then the union of
      !! every column's separable densities goes through **one**
      !! `two_electron_deriv_many` call, and every column's active two-body
      !! Gamma through **one** `active_two_electron_gradient_stacked` call. The
      !! core-Hamiltonian derivative (`build_core_hamiltonian_derivative`) and
      !! the nuclear-repulsion term are shared the same way.
      !!
      !! `cg_iterations`/`cg_residual` are indexed roots first, then extras.
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: roots(:)
      real(dp), allocatable, intent(out) :: gradients(:, :, :)   !! (3, natm, size(roots))
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: cg_tol
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations(:)   !! (size(roots) + n_extra)
      real(dp), intent(out), optional :: cg_residual(:)     !! (size(roots) + n_extra)
      real(dp), intent(in), optional :: extra_rhs(:, :)               !! (n_param, n_extra)
      real(dp), intent(in), optional :: extra_gamma(:, :, :, :, :)    !! (n_active^4, n_extra)
      real(dp), intent(in), optional :: extra_d_active(:, :, :)       !! (n_ao, n_ao, n_extra)
      real(dp), intent(in), optional :: extra_weighted(:, :, :)       !! (n_ao, n_ao, n_extra)
      real(dp), allocatable, intent(out), optional :: extra_out(:, :, :)  !! (3, natm, n_extra)

      real(dp), allocatable :: dm1_all(:, :, :), dm2_all(:, :, :, :, :)
      real(dp), allocatable :: dm1_i(:, :), dm2_i(:, :, :, :)
      real(dp), allocatable :: grad_full(:, :)
      type(mcscf_fock_t), allocatable :: fock_all(:)
      real(dp), allocatable :: rhs(:, :), x(:, :)
      real(dp), allocatable :: kappa_bar_all(:, :, :), xbar_all(:, :, :, :)
      integer, allocatable :: iterations(:)
      real(dp), allocatable :: residual(:)
      real(dp) :: use_tol
      integer :: n_states, n_roots, n_extra, n_col, use_max_iter, natm
      integer :: n_ao, n_mo, na, nb, n_param, l, j, ir, ic, ie, state_index, seg0, seg1

      ! Cheap (no derivative-integral) per-column pieces.
      real(dp), allocatable :: c_active(:, :), c_active_scratch(:, :)
      real(dp), allocatable :: cbar_active_all(:, :, :), cbar_active_scratch(:, :)
      real(dp), allocatable :: d_core(:, :), d_active_sa(:, :)
      real(dp), allocatable :: d_core_bar_all(:, :, :), d_active_bar_all(:, :, :)
      real(dp), allocatable :: weighted_bar_all(:, :, :)
      real(dp), allocatable :: d_core_bar_scratch(:, :), d_active_bar_scratch(:, :)
      real(dp), allocatable :: weighted_bar_scratch(:, :)
      real(dp), allocatable :: tdm1_scratch(:, :), tdm2_scratch(:, :, :, :)
      real(dp), allocatable :: tdm2_total_all(:, :, :, :, :)
      real(dp), allocatable :: d_ci_active_scratch(:, :), weighted_ci_scratch(:, :)
      real(dp), allocatable :: d_ci_active_all(:, :, :), weighted_ci_all(:, :, :)
      real(dp), allocatable :: ddm2_i(:, :, :, :), ddm2_all(:, :, :, :, :)
      real(dp), allocatable :: d_active_all(:, :, :), weighted_base_all(:, :, :)

      ! Shared derivative-integral quantities.
      real(dp), allocatable :: s1(:, :, :), kin(:, :, :), h1(:, :, :)
      real(dp), allocatable :: hcore_all(:, :, :, :)
      real(dp), allocatable :: density_stack(:, :, :), vhf_stack(:, :, :, :)
      real(dp), allocatable :: nuc_grad(:, :), column_grads(:, :, :)
      integer :: n_stack, slot, s_bar, s_core_bar

      if (error%has_error()) return
      n_states = size(weights)
      n_roots = size(roots)
      n_extra = 0
      if (present(extra_rhs)) n_extra = size(extra_rhs, 2)
      n_col = n_roots + n_extra
      natm = state%mol%natm
      if (present(cg_iterations)) cg_iterations = 0
      if (present(cg_residual)) cg_residual = 0.0_dp

      call check_sa_request(n_states, weights, roots, error)
      if (error%has_error()) return

      allocate (gradients(3, natm, n_roots))
      gradients = 0.0_dp
      if (present(extra_out)) then
         allocate (extra_out(3, natm, n_extra))
         extra_out = 0.0_dp
      end if
      if (n_col == 0) return

      use_tol = DEFAULT_CG_TOL
      if (present(cg_tol)) use_tol = cg_tol
      use_max_iter = DEFAULT_CG_MAX_ITER
      if (present(cg_max_iter)) use_max_iter = cg_max_iter

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)

      ! ---- one right-hand side per column: a root's own orbital gradient at
      ! the SA orbitals (its CI part is identically zero, see
      ! `sa_casscf_gradient_general`), then the caller's extra columns -------
      n_param = sa_hessian_n_param(state)
      allocate (rhs(n_param, n_col))
      rhs = 0.0_dp
      allocate (dm1_all(n_active, n_active, n_roots))
      allocate (dm2_all(n_active, n_active, n_active, n_active, n_roots))
      allocate (fock_all(n_roots))
      do ir = 1, n_roots
         state_index = roots(ir)
         ! `active_space_rdms`'s `dm1`/`dm2` are ALLOCATABLE, INTENT(OUT): an
         ! array section of another allocatable cannot itself be the actual
         ! argument, hence the per-root temporaries copied into the stack.
         call active_space_rdms(state%ci_vectors(:, :, state_index), state%alpha, state%beta, &
                                dm1_i, dm2_i, error)
         if (error%has_error()) return
         dm1_all(:, :, ir) = dm1_i
         dm2_all(:, :, :, :, ir) = dm2_i
         ! Only `general` is read from here on, and the SA Hessian state
         ! already holds the MO integrals it needs: no AO pass.
         call cheap_generalized_fock(state, dm1_i, dm2_i, fock_all(ir)%general)
         call orbital_gradient(fock_all(ir), n_inactive, n_active, grad_full)
         do l = 1, state%n_rot
            rhs(l, ir) = -grad_full(state%rows(l), state%cols(l))
         end do
      end do
      if (n_extra > 0) rhs(:, n_roots + 1:n_col) = extra_rhs

      allocate (iterations(n_col), residual(n_col))
      call sa_block_zvector_solve(state, rhs, x, iterations, residual, use_tol, &
                                  use_max_iter, error)
      if (present(cg_iterations)) cg_iterations = iterations
      if (present(cg_residual)) cg_residual = residual
      if (error%has_error()) return

      na = state%alpha%n_strings
      nb = state%beta%n_strings

      ! ---- kappa_bar/xbar for every column, from the block solution -------
      allocate (kappa_bar_all(n_mo, n_mo, n_col))
      allocate (xbar_all(na, nb, n_states, n_col))
      do ic = 1, n_col
         kappa_bar_all(:, :, ic) = 0.0_dp
         do l = 1, state%n_rot
            kappa_bar_all(state%rows(l), state%cols(l), ic) = x(l, ic)
            kappa_bar_all(state%cols(l), state%rows(l), ic) = -x(l, ic)
         end do
         do j = 1, n_states
            seg0 = state%n_rot + (j - 1)*state%n_det + 1
            seg1 = state%n_rot + j*state%n_det
            xbar_all(:, :, j, ic) = reshape(x(seg0:seg1, ic), [na, nb])
         end do
      end do

      ! ---- every column's cheap densities/weighted matrices, no derivative
      ! integral touched yet. `c_active`, `d_core` and `d_active_sa` come back
      ! identical every iteration (they depend on the SA state, not on the
      ! column), so only the last is kept. ----------------------------------
      allocate (cbar_active_all(n_ao, n_active, n_col))
      allocate (d_core_bar_all(n_ao, n_ao, n_col), d_active_bar_all(n_ao, n_ao, n_col))
      allocate (weighted_bar_all(n_ao, n_ao, n_col))
      allocate (tdm2_total_all(n_active, n_active, n_active, n_active, n_col))
      allocate (d_ci_active_all(n_ao, n_ao, n_col), weighted_ci_all(n_ao, n_ao, n_col))
      allocate (ddm2_all(n_active, n_active, n_active, n_active, n_col))
      allocate (d_active_all(n_ao, n_ao, n_roots), weighted_base_all(n_ao, n_ao, n_roots))

      do ic = 1, n_col
         call orbital_response_pieces(orbitals, n_inactive, n_active, state, &
                                      kappa_bar_all(:, :, ic), c_active_scratch, &
                                      cbar_active_scratch, d_core, d_active_sa, &
                                      d_core_bar_scratch, d_active_bar_scratch, &
                                      weighted_bar_scratch)
         cbar_active_all(:, :, ic) = cbar_active_scratch
         d_core_bar_all(:, :, ic) = d_core_bar_scratch
         d_active_bar_all(:, :, ic) = d_active_bar_scratch
         weighted_bar_all(:, :, ic) = weighted_bar_scratch

         call ci_response_pieces(orbitals, n_inactive, n_active, state, &
                                 xbar_all(:, :, :, ic), weights, c_active_scratch, d_core, &
                                 d_active_sa, tdm1_scratch, tdm2_scratch, &
                                 d_ci_active_scratch, weighted_ci_scratch, error)
         if (error%has_error()) return
         tdm2_total_all(:, :, :, :, ic) = tdm2_scratch
         d_ci_active_all(:, :, ic) = d_ci_active_scratch
         weighted_ci_all(:, :, ic) = weighted_ci_scratch

         if (ic <= n_roots) then
            call cumulant_two_particle_density(dm1_all(:, :, ic), dm2_all(:, :, :, :, ic), ddm2_i)
            ddm2_all(:, :, :, :, ic) = ddm2_i
            call build_active_density(c_active_scratch, dm1_all(:, :, ic), &
                                      d_active_all(:, :, ic))
            call build_weighted_from_fock(orbitals, fock_all(ic), weighted_base_all(:, :, ic))
         else
            ! The extra column's base term has the CI-response piece's form,
            ! and `response_separable_assemble` is linear in it, so it rides
            ! in the same slot.
            ie = ic - n_roots
            ddm2_all(:, :, :, :, ic) = extra_gamma(:, :, :, :, ie)
            d_ci_active_all(:, :, ic) = d_ci_active_all(:, :, ic) + extra_d_active(:, :, ie)
            weighted_ci_all(:, :, ic) = weighted_ci_all(:, :, ic) + extra_weighted(:, :, ie)
         end if
      end do
      c_active = c_active_scratch

      ! ---- one derivative-integral sweep for every column and separable piece
      call one_electron_deriv(state%mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(state%mol, kin, DERIV_KIN)
      call one_electron_deriv(state%mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)
      call build_core_hamiltonian_derivative(state%mol, h1, hcore_all)

      ! The orbital- and CI-response active densities meet the potentials
      ! only through `response_separable_assemble`, which is linear in them
      ! and dots their potential against `d_core` alone: one slot holds both.
      d_active_bar_all = d_active_bar_all + d_ci_active_all
      weighted_bar_all = weighted_bar_all + weighted_ci_all

      ! Slots: the two shared reference densities; three per root (its own
      ! active density, then the response's active and core parts); two per
      ! extra column (no own-density slot).
      n_stack = 2 + 3*n_roots + 2*n_extra
      allocate (density_stack(n_ao, n_ao, n_stack))
      density_stack(:, :, 1) = d_core
      density_stack(:, :, 2) = d_active_sa
      do ic = 1, n_col
         call column_slots(ic, n_roots, slot, s_bar, s_core_bar)
         if (ic <= n_roots) density_stack(:, :, slot) = d_active_all(:, :, ic)
         density_stack(:, :, s_bar) = d_active_bar_all(:, :, ic)
         density_stack(:, :, s_core_bar) = d_core_bar_all(:, :, ic)
      end do
      call two_electron_deriv_many(state%mol, density_stack, vhf_stack, error)
      deallocate (density_stack)
      if (error%has_error()) return

      allocate (nuc_grad(3, natm))
      nuc_grad = 0.0_dp
      call nuclear_repulsion_gradient(state%mol, nuc_grad)
      allocate (column_grads(3, natm, n_col))
      column_grads = 0.0_dp

      do ic = 1, n_col
         call column_slots(ic, n_roots, slot, s_bar, s_core_bar)
         if (ic <= n_roots) then
            column_grads(:, :, ic) = nuc_grad
            call base_gradient_assemble(state%mol, hcore_all, s1, vhf_stack(:, :, :, 1), &
                                        vhf_stack(:, :, :, slot), d_core, &
                                        d_active_all(:, :, ic), weighted_base_all(:, :, ic), &
                                        column_grads(:, :, ic), error)
         end if
         call response_separable_assemble(state%mol, hcore_all, s1, d_core, d_active_sa, &
                                          d_active_bar_all(:, :, ic), &
                                          weighted_bar_all(:, :, ic), vhf_stack(:, :, :, 1), &
                                          vhf_stack(:, :, :, 2), vhf_stack(:, :, :, s_bar), &
                                          column_grads(:, :, ic), error, &
                                          d_core_resp=d_core_bar_all(:, :, ic), &
                                          vhf_core_resp=vhf_stack(:, :, :, s_core_bar))
         if (error%has_error()) return
      end do

      ! ---- one stacked sweep for every column's active two-body Gamma ------
      call active_two_electron_gradient_stacked(state%mol, c_active, ddm2_all, cbar_active_all, &
                                                state%dm2_sa, tdm2_total_all, column_grads, error)
      if (error%has_error()) return

      gradients = column_grads(:, :, 1:n_roots)
      if (present(extra_out)) extra_out = column_grads(:, :, n_roots + 1:n_col)
   end subroutine sa_gradients_on_state

   pure subroutine column_slots(ic, n_roots, own, s_bar, s_core_bar)
      !! Column `ic`'s slots in `sa_gradients_on_state`'s density stack; `own`
      !! is meaningful for a root column only
      integer, intent(in) :: ic, n_roots
      integer, intent(out) :: own, s_bar, s_core_bar

      if (ic <= n_roots) then
         own = 2 + 3*(ic - 1) + 1
         s_bar = own + 1
      else
         own = 0
         s_bar = 2 + 3*n_roots + 2*(ic - n_roots - 1) + 1
      end if
      s_core_bar = s_bar + 1
   end subroutine column_slots

   subroutine response_separable_gradient(mol, d_core_ref, d_active_ref, d_active_resp, &
                                          weighted_resp, gradient, error, d_core_resp)
      !! The h1-derivative, separable two-electron-derivative and overlap
      !! (Pulay) contributions of a reference/response density pair, added
      !! into `gradient`
      !!
      !! `weighted_resp` is the caller's already-assembled, already-symmetrised
      !! AO energy-weighted matrix for this response piece (`orbital_response_gradient`'s
      !! `weighted_bar`, `ci_response_gradient`'s `weighted_ci`).
      !!
      !! `d_core_resp` absent is the CI-response case (a transition density
      !! lives on the active space alone); present, the orbital-response case.
      !! Every pairing of `{core, active} x {reference, response}` two-electron
      !! derivative-integral contraction appears **except
      !! active-response-with-active-response**, which the caller's active
      !! two-electron cumulant contraction supplies instead.
      !!
      !! Builds its own one- and two-electron derivative integrals (one
      !! `two_electron_deriv_many` sweep) and hands the per-atom assembly to
      !! `response_separable_assemble`; `sa_casscf_gradients_general` calls
      !! that one directly instead, with the sweep shared across every root.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: d_core_ref(:, :), d_active_ref(:, :)
      real(dp), intent(in) :: d_active_resp(:, :)
      real(dp), intent(in) :: weighted_resp(:, :)     !! (n_ao, n_ao)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: d_core_resp(:, :)

      real(dp), allocatable :: s1(:, :, :), kin(:, :, :), h1(:, :, :)
      real(dp), allocatable :: hcore_all(:, :, :, :)
      real(dp), allocatable :: vhf_core_ref(:, :, :), vhf_active_ref(:, :, :)
      real(dp), allocatable :: vhf_core_resp(:, :, :), vhf_active_resp(:, :, :)
      real(dp), allocatable :: density_stack(:, :, :), vhf_stack(:, :, :, :)
      integer :: n_ao, n_set
      logical :: have_core_resp

      if (error%has_error()) return
      n_ao = size(d_core_ref, 1)
      have_core_resp = present(d_core_resp)

      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)
      call build_core_hamiltonian_derivative(mol, h1, hcore_all)

      ! `d_core_ref`, `d_active_ref` and `d_active_resp` (and `d_core_resp`,
      ! present or not) each want their own `J - K/2`, with no cross term
      ! between them: one `two_electron_deriv_many` pass instead of three or
      ! four single-density ones.
      n_set = merge(4, 3, have_core_resp)
      allocate (density_stack(n_ao, n_ao, n_set))
      density_stack(:, :, 1) = d_core_ref
      density_stack(:, :, 2) = d_active_ref
      density_stack(:, :, 3) = d_active_resp
      if (have_core_resp) density_stack(:, :, 4) = d_core_resp
      call two_electron_deriv_many(mol, density_stack, vhf_stack, error)
      if (error%has_error()) return
      vhf_core_ref = vhf_stack(:, :, :, 1)
      vhf_active_ref = vhf_stack(:, :, :, 2)
      vhf_active_resp = vhf_stack(:, :, :, 3)
      if (have_core_resp) vhf_core_resp = vhf_stack(:, :, :, 4)
      deallocate (density_stack, vhf_stack)

      call response_separable_assemble(mol, hcore_all, s1, d_core_ref, d_active_ref, &
                                       d_active_resp, weighted_resp, vhf_core_ref, &
                                       vhf_active_ref, vhf_active_resp, gradient, error, &
                                       d_core_resp=d_core_resp, vhf_core_resp=vhf_core_resp)

      deallocate (s1, h1, hcore_all)
      deallocate (vhf_core_ref, vhf_active_ref, vhf_active_resp)
      if (allocated(vhf_core_resp)) deallocate (vhf_core_resp)
   end subroutine response_separable_gradient

   subroutine response_separable_assemble(mol, hcore_all, s1, d_core_ref, d_active_ref, &
                                          d_active_resp, weighted_resp, vhf_core_ref, &
                                          vhf_active_ref, vhf_active_resp, gradient, error, &
                                          d_core_resp, vhf_core_resp)
      !! `response_separable_gradient`'s per-atom assembly, given the
      !! core-Hamiltonian derivative and the two-electron potentials already
      !! built. Splitting the sweep out of the assembly is what lets
      !! `sa_casscf_gradients_general` share one `hcore_all` and one
      !! `two_electron_deriv_many` call across every root and both response
      !! pieces, instead of each root paying for three separate sweeps.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: hcore_all(:, :, :, :)   !! (n_ao, n_ao, 3, natm)
      real(dp), intent(in) :: s1(:, :, :)
      real(dp), intent(in) :: d_core_ref(:, :), d_active_ref(:, :)
      real(dp), intent(in) :: d_active_resp(:, :)
      real(dp), intent(in) :: weighted_resp(:, :)
      real(dp), intent(in) :: vhf_core_ref(:, :, :), vhf_active_ref(:, :, :)
      real(dp), intent(in) :: vhf_active_resp(:, :, :)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: d_core_resp(:, :)
      real(dp), intent(in), optional :: vhf_core_resp(:, :, :)

      real(dp), allocatable :: d_resp_total(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, iatom, comp, p0, p1
      logical :: have_core_resp

      if (error%has_error()) return
      n_ao = size(d_core_ref, 1)
      have_core_resp = present(d_core_resp)

      allocate (d_resp_total(n_ao, n_ao))
      d_resp_total = d_active_resp
      if (have_core_resp) d_resp_total = d_resp_total + d_core_resp

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(hcore_all(:, :, comp, iatom)*d_resp_total)
         end do

         if (counts(iatom) == 0) cycle

         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + 2.0_dp*sum(vhf_core_ref(p0:p1, :, comp)*d_resp_total(p0:p1, :)) &
                                    + 2.0_dp*sum(vhf_active_resp(p0:p1, :, comp)*d_core_ref(p0:p1, :)) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*weighted_resp(p0:p1, :))
            if (have_core_resp) then
               gradient(comp, iatom) = gradient(comp, iatom) &
                                       + 2.0_dp*sum(vhf_core_resp(p0:p1, :, comp)* &
                                                    (d_core_ref(p0:p1, :) + d_active_ref(p0:p1, :))) &
                                       + 2.0_dp*sum(vhf_active_ref(p0:p1, :, comp)*d_core_resp(p0:p1, :))
            end if
         end do
      end do

      deallocate (d_resp_total, offsets, counts)
   end subroutine response_separable_assemble

   subroutine build_core_hamiltonian_derivative(mol, h1, hcore_all)
      !! `response_separable_gradient`'s (and `czt_mcscf_gradient`'s) per-atom
      !! `hcore_a`, built for every atom at once. It depends only on the
      !! molecule and the density-independent core-Hamiltonian derivative
      !! `h1`, never on which density it will end up dotted against, so
      !! `sa_casscf_gradients_general` builds it exactly once and shares it
      !! across every root and every response piece.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: h1(:, :, :)
      real(dp), allocatable, intent(out) :: hcore_all(:, :, :, :)   !! (n_ao, n_ao, 3, natm)

      real(dp), allocatable :: vrinv(:, :, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, iatom, comp, p0, p1

      ! TODO(mqc): holding every atom's block costs n_ao^2 * 3 * n_atoms reals
      ! (about 300 MB at 500 AOs and 50 atoms). Looping atoms outermost in the
      ! fused assembly and contracting every density against one atom's block
      ! at a time would keep one block.
      n_ao = size(h1, 1)
      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (hcore_all(n_ao, n_ao, 3, mol%natm), vrinv(n_ao, n_ao, 3))

      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         call iprinv_deriv_at(mol, iatom, vrinv)
         vrinv = -mol%charges(iatom)*vrinv
         if (counts(iatom) > 0) then
            vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         end if
         do comp = 1, 3
            hcore_all(:, :, comp, iatom) = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
         end do
      end do

      deallocate (vrinv, offsets, counts)
   end subroutine build_core_hamiltonian_derivative

   subroutine base_gradient_assemble(mol, hcore_all, s1, vhf_core, vhf_active, d_core, &
                                     d_active, weighted, gradient, error)
      !! `czt_mcscf_gradient`'s own per-atom Hcore/Coulomb/Pulay assembly (its
      !! active two-body cumulant term is handled separately, by
      !! `active_two_electron_gradient_stacked`), from the shared
      !! `hcore_all`/`s1` and the two-electron potentials of `d_core` and this
      !! root's own active density `d_active` -- `V(d_core + d_active) =
      !! V(d_core) + V(d_active)` by linearity of `J - K/2`, so the two need
      !! not be added before contracting, only after.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: hcore_all(:, :, :, :), s1(:, :, :)
      real(dp), intent(in) :: vhf_core(:, :, :), vhf_active(:, :, :)
      real(dp), intent(in) :: d_core(:, :), d_active(:, :)
      real(dp), intent(in) :: weighted(:, :)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: d_total(:, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, iatom, comp, p0, p1

      if (error%has_error()) return
      n_ao = size(d_core, 1)
      allocate (d_total(n_ao, n_ao))
      d_total = d_core + d_active

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)

      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(hcore_all(:, :, comp, iatom)*d_total)
         end do

         if (counts(iatom) == 0) cycle

         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + 2.0_dp*sum(vhf_core(p0:p1, :, comp)*d_total(p0:p1, :)) &
                                    + 2.0_dp*sum(vhf_active(p0:p1, :, comp)*d_total(p0:p1, :)) &
                                    - 2.0_dp*sum(s1(p0:p1, :, comp)*weighted(p0:p1, :))
         end do
      end do

      deallocate (d_total, offsets, counts)
   end subroutine base_gradient_assemble

   subroutine active_two_electron_gradient_stacked(mol, c_active, ddm2_all, cbar_active_all, &
                                                   dm2_sa, tdm2_total_all, gradients, error)
      !! Every root's active two-body Gamma -- the base cumulant, the
      !! orbital-response four-leg sum, and the CI-response transition
      !! density -- summed in the AO basis one shell block at a time and
      !! contracted against the derivative integrals in **one** sweep for
      !! every root together (`active_two_electron_gradient_many`), instead of
      !! `n_root` separate sweeps (one per root, itself already one sweep
      !! since `active_two_electron_gradient_response`'s Part-A fusion).
      !! `gamma_block`'s four gemms per term are bounded by `n_active` and
      !! cost little; what this saves is the shell-quartet loop, which does
      !! not screen (`two_electron_mp2_terms`'s `with_gamma` branch has no
      !! Schwarz bound for a general four-index density) and so is paid in
      !! full every time it runs.
      !!
      !! Gamma is built only over the AOs where some leg (`c_active` or any
      !! root's `cbar_active_all`) has an amplitude above `SIGNIFICANCE` times
      !! the largest. The rest contribute exactly nothing to Gamma, and in a
      !! planar pi system that is every sigma-type AO.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: c_active(:, :)                !! (n_ao, n_active)
      real(dp), intent(in) :: ddm2_all(:, :, :, :, :)       !! (n_active^4, n_root), base cumulants
      real(dp), intent(in) :: cbar_active_all(:, :, :)      !! (n_ao, n_active, n_root)
      real(dp), intent(in) :: dm2_sa(:, :, :, :)            !! The shared SA active density
      real(dp), intent(in) :: tdm2_total_all(:, :, :, :, :)  !! (n_active^4, n_root)
      real(dp), intent(inout) :: gradients(:, :, :)         !! (3, natm, n_root), accumulated into
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: gamma_stack(:, :, :, :, :), tmp(:, :, :, :)
      real(dp), allocatable :: amplitude(:), c_sig(:, :), cbar_sig(:, :, :)
      real(dp), allocatable :: c_both(:, :), g_both(:, :, :, :)
      integer, allocatable :: ao_map(:), shells(:), sig_aos(:)
      real(dp) :: cut
      integer :: n_ao, n_root, n_sig, n_sh, sh, mu, i0, dsh, ir, n_act
      integer :: ip, ip_lo, ip_hi, p_lo, p_hi, np, per_block
      real(dp), parameter :: BLOCK_TARGET = 2.0e8_dp
      real(dp), parameter :: SIGNIFICANCE = 1.0e-12_dp

      if (error%has_error()) return
      n_ao = size(c_active, 1)
      n_root = size(ddm2_all, 5)

      allocate (amplitude(n_ao), ao_map(n_ao))
      do mu = 1, n_ao
         amplitude(mu) = max(maxval(abs(c_active(mu, :))), maxval(abs(cbar_active_all(mu, :, :))))
      end do
      cut = SIGNIFICANCE*max(maxval(amplitude), tiny(1.0_dp))
      ao_map = 0
      n_sig = 0
      do mu = 1, n_ao
         if (amplitude(mu) > cut) then
            n_sig = n_sig + 1
            ao_map(mu) = n_sig
         end if
      end do
      if (n_sig == 0) return
      allocate (sig_aos(n_sig))
      sig_aos = pack([(mu, mu=1, n_ao)], ao_map > 0)
      c_sig = c_active(sig_aos, :)
      cbar_sig = cbar_active_all(sig_aos, :, :)
      n_act = size(c_active, 2)
      allocate (c_both(n_sig, 2*n_act), g_both(2*n_act, 2*n_act, 2*n_act, 2*n_act))

      n_sh = 0
      allocate (shells(mol%nbas))
      do sh = 1, mol%nbas
         i0 = mol%shell_offset(sh)
         dsh = shell_dim(mol%cartesian, sh - 1, mol%bas)
         if (any(ao_map(i0 + 1:i0 + dsh) > 0)) then
            n_sh = n_sh + 1
            shells(n_sh) = sh
         end if
      end do
      shells = shells(1:n_sh)

      per_block = max(1, int(BLOCK_TARGET/(2.0_dp*real(n_sig, dp)**3*8.0_dp*real(n_root, dp))))

      ip_lo = 1
      do while (ip_lo <= n_sh)
         p_lo = first_mapped(shells(ip_lo))
         ip_hi = ip_lo
         do ip = ip_lo, n_sh
            if (ip > ip_lo .and. last_mapped(shells(ip)) - p_lo + 1 > per_block) exit
            ip_hi = ip
         end do
         p_hi = last_mapped(shells(ip_hi))
         np = p_hi - p_lo + 1

         allocate (gamma_stack(np, n_sig, n_sig, n_sig, n_root))
         do ir = 1, n_root
            ! One transform over the doubled leg `[C | Cbar]` gives all six
            ! terms: the base and CI-response densities in the all-`C` block,
            ! `dm2_sa` in each block with `Cbar` on exactly one leg.
            c_both(:, 1:n_act) = c_sig
            c_both(:, n_act + 1:2*n_act) = cbar_sig(:, :, ir)
            g_both = 0.0_dp
            g_both(1:n_act, 1:n_act, 1:n_act, 1:n_act) = ddm2_all(:, :, :, :, ir) &
                                                         + tdm2_total_all(:, :, :, :, ir)
            g_both(n_act + 1:, 1:n_act, 1:n_act, 1:n_act) = dm2_sa
            g_both(1:n_act, n_act + 1:, 1:n_act, 1:n_act) = dm2_sa
            g_both(1:n_act, 1:n_act, n_act + 1:, 1:n_act) = dm2_sa
            g_both(1:n_act, 1:n_act, 1:n_act, n_act + 1:) = dm2_sa
            call gamma_block(c_both, g_both, p_lo, p_hi, tmp)
            gamma_stack(:, :, :, :, ir) = tmp
         end do

         call active_two_electron_gradient_many(mol, gamma_stack, shells, ip_lo, ip_hi, ao_map, &
                                                p_lo - 1, gradients, error)
         deallocate (gamma_stack)
         if (error%has_error()) return
         ip_lo = ip_hi + 1
      end do

   contains

      function first_mapped(shell) result(first)
         !! The smallest compressed index among `shell`'s AOs
         integer, intent(in) :: shell
         integer :: first
         integer :: q
         first = huge(1)
         do q = mol%shell_offset(shell) + 1, mol%shell_offset(shell) + &
            shell_dim(mol%cartesian, shell - 1, mol%bas)
            if (ao_map(q) > 0) first = min(first, ao_map(q))
         end do
      end function first_mapped

      function last_mapped(shell) result(last)
         !! The largest compressed index among `shell`'s AOs
         integer, intent(in) :: shell
         integer :: last
         integer :: q
         last = 0
         do q = mol%shell_offset(shell) + 1, mol%shell_offset(shell) + &
            shell_dim(mol%cartesian, shell - 1, mol%bas)
            last = max(last, ao_map(q))
         end do
      end function last_mapped
   end subroutine active_two_electron_gradient_stacked

   subroutine build_active_density(c_active, dm1, d_active)
      !! The AO active density `c_active dm1 c_active^T` -- the base
      !! gradient's own `d_active`, built once per root here from that root's
      !! own `dm1` rather than the SA one `orbital_response_pieces` builds.
      real(dp), intent(in) :: c_active(:, :)     !! (n_ao, n_active)
      real(dp), intent(in) :: dm1(:, :)          !! (n_active, n_active)
      real(dp), intent(out) :: d_active(:, :)    !! (n_ao, n_ao)

      real(dp), allocatable :: work(:, :)

      allocate (work(size(c_active, 1), size(c_active, 2)))
      call pic_gemm(c_active, dm1, work, beta=0.0_dp)
      call pic_gemm(work, c_active, d_active, transb="T", beta=0.0_dp)
      deallocate (work)
   end subroutine build_active_density

   subroutine build_weighted_from_fock(orbitals, fock, weighted)
      !! `czt_mcscf_gradient`'s own energy-weighted density,
      !! `C sym(F) C^T`, from a root's own generalised Fock -- built once per
      !! root here so the base term's Pulay contribution can be assembled
      !! from a shared, precomputed overlap derivative later.
      real(dp), intent(in) :: orbitals(:, :)
      type(mcscf_fock_t), intent(in) :: fock
      real(dp), intent(out) :: weighted(:, :)   !! (n_ao, n_ao)

      real(dp), allocatable :: work(:, :)

      allocate (work(size(orbitals, 1), size(orbitals, 2)))
      call pic_gemm(orbitals, 0.5_dp*(fock%general + transpose(fock%general)), work, &
                    beta=0.0_dp)
      call pic_gemm(work, orbitals, weighted, transb="T", beta=0.0_dp)
      deallocate (work)
   end subroutine build_weighted_from_fock

   subroutine orbital_response_gradient(mol, orbitals, n_inactive, n_active, state, &
                                        kappa_bar, gradient, error)
      !! `d/dR[kappa_bar . dE_SA/dkappa]`, at fixed `(dm1_sa, dm2_sa)` -- see
      !! the module docstring for the term-by-term derivation
      !!
      !! `orbital_response_pieces` builds every density and energy-weighted
      !! matrix this needs (no derivative integral touched yet); the two
      !! calls after it are the only derivative-integral sweeps this term
      !! costs. The cross-root fused path (`sa_casscf_gradients_general`)
      !! calls the same builder per root and shares those two sweeps across
      !! every root instead of repeating them here.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: kappa_bar(:, :)     !! (n_mo, n_mo), antisymmetric
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: c_active(:, :), cbar_active(:, :)
      real(dp), allocatable :: d_core(:, :), d_active(:, :)
      real(dp), allocatable :: d_core_bar(:, :), d_active_bar(:, :)
      real(dp), allocatable :: weighted_bar(:, :)
      real(dp), allocatable :: active_2body(:, :)

      if (error%has_error()) return
      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp

      call orbital_response_pieces(orbitals, n_inactive, n_active, state, kappa_bar, &
                                   c_active, cbar_active, d_core, d_active, d_core_bar, &
                                   d_active_bar, weighted_bar)

      call response_separable_gradient(mol, d_core, d_active, d_active_bar, weighted_bar, &
                                       gradient, error, d_core_resp=d_core_bar)
      if (error%has_error()) return

      ! The active two-body response: the four-term one-index transform of
      ! the FULL active two-particle density `dm2_sa`, summed over all four
      ! legs before the derivative-integral sweep rather than one leg per
      ! sweep (`active_two_electron_gradient_response`). Unlike the base
      ! gradient's own *energy* -- which splits the active-active
      ! two-electron contribution into a separable "classical" mean-field
      ! piece (folded into `response_separable_gradient` above via `d_total`)
      ! plus a cumulant correction (`ddm2`) -- the *generalized Fock* this
      ! routine differentiates never makes that split: `generalized_fock`'s
      ! active-row block contracts the full `dm2_sa` against `eri_gaaa`
      ! directly (its inactive-row block is where a classical active mean
      ! field appears, already covered by `d_active`/`d_active_bar` above).
      ! Using the cumulant here instead of the full density under-counts this
      ! term by roughly an order of magnitude on LiH/STO-3G SA-2 -- caught by
      ! `test_orbital_response_connection` together with an isolated
      ! fixed-orbital/dm2-zeroed bisection that pinned the discrepancy to
      ! exactly this piece.
      allocate (active_2body(3, mol%natm))
      active_2body = 0.0_dp
      call active_two_electron_gradient_response(mol, c_active, cbar_active, state%dm2_sa, &
                                                 active_2body, error)
      if (error%has_error()) return
      gradient = gradient + active_2body
   end subroutine orbital_response_gradient

   subroutine orbital_response_pieces(orbitals, n_inactive, n_active, state, kappa_bar, &
                                      c_active, cbar_active, d_core, d_active, d_core_bar, &
                                      d_active_bar, weighted_bar)
      !! Every density and energy-weighted matrix `orbital_response_gradient`
      !! needs, built with no derivative integral: the one-index-transformed
      !! response densities from `kappa_bar`, and the overlap term's `W`
      !! matrix. Shared by the single-root path and, per root, by
      !! `sa_casscf_gradients_general`'s cross-root fusion.
      !!
      !! **The overlap/Pulay term.** Exactly `czt_mcscf_gradient`'s own
      !! overlap term, generalised from a stationary energy to this
      !! non-stationary Lagrange term: `kappa_bar . dE_SA/dkappa` depends on
      !! `C` both explicitly (through the AO integrals, at fixed
      !! `dm1_sa`/`dm2_sa`, handled by the caller's derivative-integral
      !! sweeps) and through the connection `C(R) = C0 (C0^T S(R) C0)^(-1/2)`
      !! that keeps a fixed numerical `C0` orthonormal as `R` moves. `W` below
      !! is the derivative of the Lagrange term wrt a *general* (not just
      !! antisymmetric) orbital change `C -> C(1+A)`, and the connection's own
      !! derivative is `C0 -> C0(1 - 1/2 dS_MO/dR)`, a *symmetric* `A`; the
      !! overlap term is `-1/2 sum_pq (dS_MO/dR)_pq (W_pq + W_qp)`, i.e.
      !! `-s1 . sym(W)`, the same shape `czt_mcscf_gradient`'s Pulay term takes
      !! with `F -> W`.
      !!
      !! `W` itself: write `H(A, eps) = E_SA(C(1+A) exp(eps kappa_bar))`. `W =
      !! d^2H/dA deps` at `(0,0)`, and Clairaut's theorem lets the two
      !! derivatives be taken in either order. Taking `d/dA` first, at fixed
      !! `eps`, needs `C(1+A) exp(eps kappa_bar) = [C exp(eps kappa_bar)] (1 +
      !! D^-1 A D)` with `D = exp(eps kappa_bar)` (inserting `D^-1 D = 1`), so
      !! `dH/dA|_{A=0} = D F_sa(C exp(eps kappa_bar)) D^T` (`D` orthogonal,
      !! `D^-1 = D^T`) by the generalised Fock's own definition as `dE/dA` at
      !! `C exp(eps kappa_bar)`. Differentiating that wrt `eps` at `eps=0`
      !! (`D=1`, `dD/deps=kappa_bar`) gives `W = F_bar + kappa_bar . F_sa -
      !! F_sa . kappa_bar`, built below as `sym_f_bar` (from `one_index_fock`)
      !! plus the commutator sandwiched between `cbar`/`orbitals` and
      !! `orbitals`/`orbitals`. Verified to `~1e-11` against
      !! `test_orbital_response_connection`'s connection-based reference, no
      !! extra factor.
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: kappa_bar(:, :)
      real(dp), allocatable, intent(out) :: c_active(:, :), cbar_active(:, :)
      real(dp), allocatable, intent(out) :: d_core(:, :), d_active(:, :)
      real(dp), allocatable, intent(out) :: d_core_bar(:, :), d_active_bar(:, :)
      real(dp), allocatable, intent(out) :: weighted_bar(:, :)

      real(dp), allocatable :: c_core(:, :)
      real(dp), allocatable :: cbar(:, :), cbar_core(:, :)
      real(dp), allocatable :: f_bar(:, :), sym_f_sa(:, :), sym_f_bar(:, :)
      real(dp), allocatable :: t1(:, :), work(:, :)
      integer :: n_ao, n_mo, n_occ

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_occ = n_inactive + n_active

      allocate (c_active(n_ao, n_active))
      c_active = orbitals(:, n_inactive + 1:n_occ)

      allocate (cbar(n_ao, n_mo))
      call pic_gemm(orbitals, kappa_bar, cbar, beta=0.0_dp)
      allocate (cbar_active(n_ao, n_active))
      cbar_active = cbar(:, n_inactive + 1:n_occ)

      allocate (d_core(n_ao, n_ao), d_active(n_ao, n_ao))
      d_core = 0.0_dp
      if (n_inactive > 0) then
         allocate (c_core(n_ao, n_inactive))
         c_core = orbitals(:, 1:n_inactive)
         call pic_gemm(c_core, c_core, d_core, transb="T", alpha=2.0_dp, beta=0.0_dp)
      end if
      d_active = 0.0_dp
      if (n_active > 0) then
         allocate (work(n_ao, n_active))
         call pic_gemm(c_active, state%dm1_sa, work, beta=0.0_dp)
         call pic_gemm(work, c_active, d_active, transb="T", beta=0.0_dp)
         deallocate (work)
      end if

      allocate (d_core_bar(n_ao, n_ao))
      d_core_bar = 0.0_dp
      if (n_inactive > 0) then
         allocate (cbar_core(n_ao, n_inactive))
         cbar_core = cbar(:, 1:n_inactive)
         call pic_gemm(cbar_core, c_core, d_core_bar, transb="T", alpha=2.0_dp, beta=0.0_dp)
         d_core_bar = d_core_bar + transpose(d_core_bar)
      end if

      allocate (d_active_bar(n_ao, n_ao))
      d_active_bar = 0.0_dp
      if (n_active > 0) then
         allocate (work(n_ao, n_active))
         call pic_gemm(cbar_active, state%dm1_sa, work, beta=0.0_dp)
         call pic_gemm(work, c_active, d_active_bar, transb="T", beta=0.0_dp)
         deallocate (work)
         d_active_bar = d_active_bar + transpose(d_active_bar)
      end if

      allocate (f_bar(n_mo, n_mo))
      call one_index_fock(state%a_block, state%b_block, state%fock_sa, state%dm1_sa, &
                          state%dm2_sa, n_inactive, n_active, kappa_bar, f_bar)
      allocate (sym_f_sa(n_mo, n_mo), sym_f_bar(n_mo, n_mo))
      sym_f_sa = 0.5_dp*(state%fock_sa%general + transpose(state%fock_sa%general))
      sym_f_bar = 0.5_dp*(f_bar + transpose(f_bar))

      allocate (t1(n_ao, n_ao), work(n_ao, n_mo))
      call pic_gemm(cbar, sym_f_sa, work, beta=0.0_dp)
      call pic_gemm(work, orbitals, t1, transb="T", beta=0.0_dp)
      allocate (weighted_bar(n_ao, n_ao))
      weighted_bar = t1 + transpose(t1)
      call pic_gemm(orbitals, sym_f_bar, work, beta=0.0_dp)
      call pic_gemm(work, orbitals, t1, transb="T", beta=0.0_dp)
      weighted_bar = weighted_bar + t1
      deallocate (t1, work)
   end subroutine orbital_response_pieces

   subroutine ci_response_gradient(mol, orbitals, n_inactive, n_active, state, xbar, &
                                   weights, gradient, error)
      !! `d/dR[sum_J xbar_J . 2 w_J (H - E_J) c_J]`, from the symmetrised
      !! transition density between `xbar(:,:,J)` and `c_J`, summed over `J`
      !!
      !! `xbar` is projected off every averaged state here, so any component
      !! the Z-vector solve leaves along them has no effect. Those components
      !! are null directions of the projected Hessian, and their Lagrange
      !! multipliers are zero for non-degenerate states.
      !!
      !! `ci_response_pieces` builds the transition densities and the
      !! energy-weighted matrix (no derivative integral touched yet); the two
      !! calls after it are the only derivative-integral sweeps this term
      !! costs, shared across roots by `sa_casscf_gradients_general`'s fused
      !! path the same way `orbital_response_gradient`'s are.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: xbar(:, :, :)      !! (na, nb, n_states)
      real(dp), intent(in) :: weights(:)
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: c_active(:, :), d_core(:, :), d_active(:, :)
      real(dp), allocatable :: tdm1_total(:, :), tdm2_total(:, :, :, :)
      real(dp), allocatable :: d_ci_active(:, :), weighted_ci(:, :)

      if (error%has_error()) return
      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp
      if (n_active == 0) return   ! no active space, nothing a CI response touches

      call ci_response_pieces(orbitals, n_inactive, n_active, state, xbar, weights, &
                              c_active, d_core, d_active, tdm1_total, tdm2_total, &
                              d_ci_active, weighted_ci, error)
      if (error%has_error()) return

      call response_separable_gradient(mol, d_core, d_active, d_ci_active, weighted_ci, &
                                       gradient, error)
      if (error%has_error()) return

      call active_two_electron_gradient(mol, c_active, tdm2_total, gradient, error)
      if (error%has_error()) return
   end subroutine ci_response_gradient

   subroutine ci_response_pieces(orbitals, n_inactive, n_active, state, xbar, weights, &
                                 c_active, d_core, d_active, tdm1_total, tdm2_total, &
                                 d_ci_active, weighted_ci, error)
      !! Every density and energy-weighted matrix `ci_response_gradient`
      !! needs, built with no derivative integral: the weighted-sum
      !! transition densities from `xbar` and every averaged state, and the
      !! energy-weighted matrix `cheap_generalized_fock` builds from them.
      !! `n_active > 0` is the caller's responsibility (`ci_response_gradient`
      !! itself returns early otherwise).
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: xbar(:, :, :)
      real(dp), intent(in) :: weights(:)
      real(dp), allocatable, intent(out) :: c_active(:, :), d_core(:, :), d_active(:, :)
      real(dp), allocatable, intent(out) :: tdm1_total(:, :), tdm2_total(:, :, :, :)
      real(dp), allocatable, intent(out) :: d_ci_active(:, :), weighted_ci(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: c_core(:, :), work(:, :)
      real(dp), allocatable :: tdm1a(:, :), tdm2a(:, :, :, :)
      real(dp), allocatable :: tdm1_j(:, :), tdm2_j(:, :, :, :)
      real(dp), allocatable :: fock_ci(:, :), sym_fock_ci(:, :)
      real(dp), allocatable :: xbar_j(:, :)
      integer :: n_ao, n_mo, n_occ, j

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_occ = n_inactive + n_active

      allocate (c_active(n_ao, n_active))
      c_active = orbitals(:, n_inactive + 1:n_occ)
      allocate (d_core(n_ao, n_ao))
      d_core = 0.0_dp
      if (n_inactive > 0) then
         allocate (c_core(n_ao, n_inactive))
         c_core = orbitals(:, 1:n_inactive)
         call pic_gemm(c_core, c_core, d_core, transb="T", alpha=2.0_dp, beta=0.0_dp)
      end if
      allocate (d_active(n_ao, n_ao))
      d_active = 0.0_dp
      allocate (work(n_ao, n_active))
      call pic_gemm(c_active, state%dm1_sa, work, beta=0.0_dp)
      call pic_gemm(work, c_active, d_active, transb="T", beta=0.0_dp)
      deallocate (work)

      do j = 1, size(weights)
         xbar_j = xbar(:, :, j)
         call project_ci_block(state, xbar_j)
         call transition_rdms(state%ci_vectors(:, :, j), xbar_j, state%alpha, &
                              state%beta, tdm1a, tdm2a, error)
         if (error%has_error()) return
         tdm1_j = tdm1a + transpose(tdm1a)
         tdm2_j = tdm2a + reshape(tdm2a, shape(tdm2a), order=[2, 1, 4, 3])
         if (j == 1) then
            tdm1_total = weights(1)*tdm1_j
            tdm2_total = weights(1)*tdm2_j
         else
            tdm1_total = tdm1_total + weights(j)*tdm1_j
            tdm2_total = tdm2_total + weights(j)*tdm2_j
         end if
      end do

      allocate (d_ci_active(n_ao, n_ao), work(n_ao, n_active))
      call pic_gemm(c_active, tdm1_total, work, beta=0.0_dp)
      call pic_gemm(work, c_active, d_ci_active, transb="T", beta=0.0_dp)
      deallocate (work)

      call cheap_generalized_fock(state, tdm1_total, tdm2_total, fock_ci, delta_only=.true.)
      allocate (sym_fock_ci(n_mo, n_mo))
      sym_fock_ci = 0.5_dp*(fock_ci + transpose(fock_ci))
      allocate (weighted_ci(n_ao, n_ao), work(n_ao, n_mo))
      call pic_gemm(orbitals, sym_fock_ci, work, beta=0.0_dp)
      call pic_gemm(work, orbitals, weighted_ci, transb="T", beta=0.0_dp)
      deallocate (work)
   end subroutine ci_response_pieces

end module mqc_czt_sa_gradient
