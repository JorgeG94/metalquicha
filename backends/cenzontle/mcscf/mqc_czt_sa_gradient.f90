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
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks
   use mqc_czt_gradient, only: one_electron_deriv, iprinv_deriv_at, two_electron_deriv, &
                               DERIV_OVLP, DERIV_KIN, DERIV_NUC
   use mqc_czt_mcscf, only: mcscf_fock_t, generalized_fock, orbital_gradient, one_index_fock
   use mqc_czt_mcscf_gradient, only: czt_mcscf_gradient, active_two_electron_gradient
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
      !! `czt_mcscf_gradient`, bit for bit, and then `state_index` must be 1.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)         !! (n_ao, n_mo), the SA orbitals
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)    !! (na, nb, >= n_states)
      real(dp), intent(in) :: energies(:)            !! (>= n_states), total
      real(dp), intent(in) :: weights(:)              !! (n_states)
      integer, intent(in) :: state_index              !! Which root, 1-based
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, n_atoms), Hartree/Bohr
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
         ! The general path's range check, repeated because this branch never
         ! reaches it: without it, any `state_index` returns root 1.
         if (state_index /= 1) then
            call error%set(ERROR_VALIDATION, "sa_casscf_gradient: state "// &
                           to_char(state_index)//" is not one of the 1 averaged "// &
                           "states.")
            return
         end if
         call build_link_table(n_active, n_alpha, alpha, error)
         if (error%has_error()) return
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
      !! `czt_sa_casscf_gradient` without its `n_states = 1` short-circuit
      !!
      !! At `n_states = 1` this reproduces `czt_mcscf_gradient` to the Z-vector
      !! solve's tolerance, not bit for bit. Arguments as there.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)         !! (n_ao, n_mo), the SA orbitals
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)    !! (na, nb, >= n_states)
      real(dp), intent(in) :: energies(:)            !! (>= n_states), total
      real(dp), intent(in) :: weights(:)             !! (n_states)
      integer, intent(in) :: state_index             !! Which root, 1-based
      real(dp), allocatable, intent(out) :: gradient(:, :)   !! (3, n_atoms), Hartree/Bohr
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
      !! Preconditioned CG on `H_SA x = rhs`, block-shaped with `n_vec = 1`
      !! columns so a later fused multi-root solve (phase 5) can widen `rhs`
      !! without changing this routine
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: rhs(:, :)      !! (n_param, n_vec)
      real(dp), allocatable, intent(out) :: x(:, :)
      integer, intent(out) :: iterations
      real(dp), intent(out) :: residual   !! The relative residual reached
      real(dp), intent(in) :: tol
      integer, intent(in) :: max_iter
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: r(:, :), z(:, :), p(:, :), ap(:, :)
      real(dp) :: rz, rz_new, pap, step, target_norm, r_norm
      integer :: n_param, n_vec, iter

      if (error%has_error()) return
      n_param = size(rhs, 1)
      n_vec = size(rhs, 2)
      allocate (x(n_param, n_vec), r(n_param, n_vec), z(n_param, n_vec))
      allocate (p(n_param, n_vec), ap(n_param, n_vec))

      x = 0.0_dp
      r = rhs
      target_norm = sqrt(sum(rhs*rhs))
      iterations = 0
      residual = 0.0_dp
      if (target_norm <= 0.0_dp) return

      call sa_hessian_precondition(state, r, z)
      p = z
      rz = sum(r*z)

      do iter = 1, max_iter
         iterations = iter
         call sa_hessian_apply(state, p, ap, error)
         if (error%has_error()) return
         pap = sum(p*ap)
         if (pap <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "sa_casscf_gradient: the Z-vector CG search "// &
                           "direction has non-positive curvature, so the SA Hessian is "// &
                           "not positive definite on the non-redundant space -- the "// &
                           "reference is not a genuine SA-CASSCF minimum.")
            return
         end if
         step = rz/pap
         x = x + step*p
         r = r - step*ap
         r_norm = sqrt(sum(r*r))
         residual = r_norm/target_norm
         if (residual <= tol) exit
         call sa_hessian_precondition(state, r, z)
         rz_new = sum(r*z)
         p = z + (rz_new/rz)*p
         rz = rz_new
      end do

      deallocate (r, z, p, ap)
   end subroutine sa_zvector_solve

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
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: d_core_ref(:, :), d_active_ref(:, :)
      real(dp), intent(in) :: d_active_resp(:, :)
      real(dp), intent(in) :: weighted_resp(:, :)     !! (n_ao, n_ao)
      real(dp), intent(inout) :: gradient(:, :)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: d_core_resp(:, :)

      real(dp), allocatable :: s1(:, :, :), kin(:, :, :), h1(:, :, :)
      real(dp), allocatable :: vrinv(:, :, :), hcore_a(:, :, :)
      real(dp), allocatable :: d_resp_total(:, :)
      real(dp), allocatable :: vhf_core_ref(:, :, :), vhf_active_ref(:, :, :)
      real(dp), allocatable :: vhf_core_resp(:, :, :), vhf_active_resp(:, :, :)
      integer, allocatable :: offsets(:), counts(:)
      integer :: n_ao, iatom, comp, p0, p1
      logical :: have_core_resp

      if (error%has_error()) return
      n_ao = size(d_core_ref, 1)
      have_core_resp = present(d_core_resp)

      allocate (d_resp_total(n_ao, n_ao))
      d_resp_total = d_active_resp
      if (have_core_resp) d_resp_total = d_resp_total + d_core_resp

      call one_electron_deriv(mol, s1, DERIV_OVLP)
      s1 = -s1
      call one_electron_deriv(mol, kin, DERIV_KIN)
      call one_electron_deriv(mol, h1, DERIV_NUC)
      h1 = -(kin + h1)
      deallocate (kin)

      call two_electron_deriv(mol, d_core_ref, vhf_core_ref, error)
      if (error%has_error()) return
      call two_electron_deriv(mol, d_active_ref, vhf_active_ref, error)
      if (error%has_error()) return
      call two_electron_deriv(mol, d_active_resp, vhf_active_resp, error)
      if (error%has_error()) return
      if (have_core_resp) then
         call two_electron_deriv(mol, d_core_resp, vhf_core_resp, error)
         if (error%has_error()) return
      end if

      allocate (offsets(mol%natm), counts(mol%natm))
      call atom_ao_blocks(mol, offsets, counts)
      allocate (vrinv(n_ao, n_ao, 3), hcore_a(n_ao, n_ao, 3))

      do iatom = 1, mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)

         call iprinv_deriv_at(mol, iatom, vrinv)
         vrinv = -mol%charges(iatom)*vrinv
         if (counts(iatom) > 0) then
            vrinv(p0:p1, :, :) = vrinv(p0:p1, :, :) + h1(p0:p1, :, :)
         end if
         do comp = 1, 3
            hcore_a(:, :, comp) = vrinv(:, :, comp) + transpose(vrinv(:, :, comp))
         end do

         do comp = 1, 3
            gradient(comp, iatom) = gradient(comp, iatom) &
                                    + sum(hcore_a(:, :, comp)*d_resp_total)
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

      deallocate (s1, h1, vrinv, hcore_a, d_resp_total, offsets, counts)
      deallocate (vhf_core_ref, vhf_active_ref, vhf_active_resp)
      if (allocated(vhf_core_resp)) deallocate (vhf_core_resp)
   end subroutine response_separable_gradient

   subroutine orbital_response_gradient(mol, orbitals, n_inactive, n_active, state, &
                                        kappa_bar, gradient, error)
      !! `d/dR[kappa_bar . dE_SA/dkappa]`, at fixed `(dm1_sa, dm2_sa)` -- see
      !! the module docstring for the term-by-term derivation
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: kappa_bar(:, :)     !! (n_mo, n_mo), antisymmetric
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: c_core(:, :), c_active(:, :)
      real(dp), allocatable :: cbar(:, :), cbar_core(:, :), cbar_active(:, :)
      real(dp), allocatable :: d_core(:, :), d_active(:, :)
      real(dp), allocatable :: d_core_bar(:, :), d_active_bar(:, :)
      real(dp), allocatable :: f_bar(:, :), sym_f_sa(:, :), sym_f_bar(:, :)
      real(dp), allocatable :: weighted_bar(:, :), t1(:, :), work(:, :)
      real(dp), allocatable :: active_2body(:, :)
      integer :: n_ao, n_mo, n_occ

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_occ = n_inactive + n_active

      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp

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

      ! ---- the overlap/Pulay term -------------------------------------------
      ! Exactly `czt_mcscf_gradient`'s own overlap term, generalised from a
      ! stationary energy to this non-stationary Lagrange term: `kappa_bar .
      ! dE_SA/dkappa` depends on `C` both explicitly (through the AO
      ! integrals, at fixed `dm1_sa`/`dm2_sa`, handled by the rest of this
      ! routine) and through the connection `C(R) = C0 (C0^T S(R) C0)^(-1/2)`
      ! that keeps a fixed numerical `C0` orthonormal as `R` moves. `W` below
      ! is the derivative of the Lagrange term wrt a *general* (not just
      ! antisymmetric) orbital change `C -> C(1+A)`, and the connection's own
      ! derivative is `C0 -> C0(1 - 1/2 dS_MO/dR)`, a *symmetric* `A`; the
      ! overlap term is `-1/2 sum_pq (dS_MO/dR)_pq (W_pq + W_qp)`, i.e.
      ! `-s1 . sym(W)`, the same shape `czt_mcscf_gradient`'s Pulay term takes
      ! with `F -> W`.
      !
      ! `W` itself: write `H(A, eps) = E_SA(C(1+A) exp(eps kappa_bar))`. `W =
      ! d^2H/dA deps` at `(0,0)`, and Clairaut's theorem lets the two
      ! derivatives be taken in either order. Taking `d/dA` first, at fixed
      ! `eps`, needs `C(1+A) exp(eps kappa_bar) = [C exp(eps kappa_bar)] (1 +
      ! D^-1 A D)` with `D = exp(eps kappa_bar)` (inserting `D^-1 D = 1`), so
      ! `dH/dA|_{A=0} = D F_sa(C exp(eps kappa_bar)) D^T` (`D` orthogonal,
      ! `D^-1 = D^T`) by the generalised Fock's own definition as `dE/dA` at
      ! `C exp(eps kappa_bar)`. Differentiating that wrt `eps` at `eps=0`
      ! (`D=1`, `dD/deps=kappa_bar`) gives `W = F_bar + kappa_bar . F_sa -
      ! F_sa . kappa_bar`, built below as `sym_f_bar` (from `one_index_fock`)
      ! plus the commutator sandwiched between `cbar`/`orbitals` and
      ! `orbitals`/`orbitals`. Verified to `~1e-11` against
      ! `test_orbital_response_connection`'s connection-based reference, no
      ! extra factor.
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

      call response_separable_gradient(mol, d_core, d_active, d_active_bar, weighted_bar, &
                                       gradient, error, d_core_resp=d_core_bar)
      if (error%has_error()) return

      ! ---- the active two-body response: the four-term one-index transform
      ! of the FULL active two-particle density `dm2_sa`, one AO-transform
      ! leg at a time. Unlike the base gradient's own *energy* -- which
      ! splits the active-active two-electron contribution into a separable
      ! "classical" mean-field piece (folded into `response_separable_gradient`
      ! above via `d_total`) plus a cumulant correction (`ddm2`) -- the
      ! *generalized Fock* this routine differentiates never makes that
      ! split: `generalized_fock`'s active-row block contracts the full
      ! `dm2_sa` against `eri_gaaa` directly (its inactive-row block is where
      ! a classical active mean field appears, already covered by `d_active`/
      ! `d_active_bar` above). Using the cumulant here instead of the full
      ! density under-counts this term by roughly an order of magnitude on
      ! LiH/STO-3G SA-2 -- caught by `test_orbital_response_connection`
      ! together with an isolated fixed-orbital/dm2-zeroed bisection that
      ! pinned the discrepancy to exactly this piece.
      allocate (active_2body(3, mol%natm))
      active_2body = 0.0_dp
      call active_two_electron_gradient(mol, cbar_active, state%dm2_sa, active_2body, error, &
                                        c2=c_active, c3=c_active, c4=c_active)
      if (error%has_error()) return
      call active_two_electron_gradient(mol, c_active, state%dm2_sa, active_2body, error, &
                                        c2=cbar_active, c3=c_active, c4=c_active)
      if (error%has_error()) return
      call active_two_electron_gradient(mol, c_active, state%dm2_sa, active_2body, error, &
                                        c2=c_active, c3=cbar_active, c4=c_active)
      if (error%has_error()) return
      call active_two_electron_gradient(mol, c_active, state%dm2_sa, active_2body, error, &
                                        c2=c_active, c3=c_active, c4=cbar_active)
      if (error%has_error()) return
      gradient = gradient + active_2body
   end subroutine orbital_response_gradient

   subroutine ci_response_gradient(mol, orbitals, n_inactive, n_active, state, xbar, &
                                   weights, gradient, error)
      !! `d/dR[sum_J xbar_J . 2 w_J (H - E_J) c_J]`, from the symmetrised
      !! transition density between `xbar(:,:,J)` and `c_J`, summed over `J`
      !!
      !! `xbar` is projected off every averaged state here, so any component
      !! the Z-vector solve leaves along them has no effect. Those components
      !! are null directions of the projected Hessian, and their Lagrange
      !! multipliers are zero for non-degenerate states.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: xbar(:, :, :)      !! (na, nb, n_states)
      real(dp), intent(in) :: weights(:)
      real(dp), allocatable, intent(out) :: gradient(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: c_core(:, :), c_active(:, :)
      real(dp), allocatable :: d_core(:, :), d_active(:, :)
      real(dp), allocatable :: tdm1a(:, :), tdm2a(:, :, :, :)
      real(dp), allocatable :: tdm1_j(:, :), tdm2_j(:, :, :, :)
      real(dp), allocatable :: tdm1_total(:, :), tdm2_total(:, :, :, :)
      real(dp), allocatable :: d_ci_active(:, :), work(:, :)
      real(dp), allocatable :: fock_ci(:, :), sym_fock_ci(:, :), weighted_ci(:, :)
      real(dp), allocatable :: xbar_j(:, :)
      integer :: n_ao, n_mo, n_occ, j

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_occ = n_inactive + n_active

      allocate (gradient(3, mol%natm))
      gradient = 0.0_dp
      if (n_active == 0) return   ! no active space, nothing a CI response touches

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

      call response_separable_gradient(mol, d_core, d_active, d_ci_active, weighted_ci, &
                                       gradient, error)
      if (error%has_error()) return

      call active_two_electron_gradient(mol, c_active, tdm2_total, gradient, error)
      if (error%has_error()) return
   end subroutine ci_response_gradient

end module mqc_czt_sa_gradient
