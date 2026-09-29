!! The matrix-free Hessian-vector product of the SA-CASSCF energy
module mqc_czt_sa_hessian
   !! Equations in `mqc_docs/source/developer_sa_casscf.rst`, "The SA-CASSCF
   !! Lagrangian and Z-vector" and "The SA Hessian blocks".
   !!
   !! **The joint parameter vector.** A point in the space the SA Hessian acts
   !! on is `[kappa (n_rot) ; x_1 (n_det) ; ... ; x_N (n_det)]`: `kappa` one
   !! entry per `rotation_parameters` pair (the orbital rotation, same order
   !! `orbital_gradient`/`orbital_hessian` use), and each `x_J` a full
   !! determinant-space vector -- the first-order change of state `J`'s CI
   !! vector, `c_J -> c_J + x_J`, with `x_J` orthogonal to `c_J` (PySCF's
   !! `newton_casscf` convention: the redundant "renormalise `c_J`" direction
   !! is not part of the parameter space at all). `sa_hessian_n_param` gives
   !! the flat length; `sa_hessian_apply`/`sa_hessian_precondition` work on
   !! that flat layout directly, block-shaped over several vectors at once.
   !!
   !! **Every quantity built once, in `build_sa_hessian`, off a converged
   !! `run_czt_casscf` result.** `a_block`/`b_block` (the `(n_mo, n_occ, n_mo,
   !! n_occ)` MO integrals `orbital_hessian` builds once per macro-iteration),
   !! the SA generalised Fock, the reference CI vectors and their active
   !! energies, and the folded active-space Hamiltonian for the CI sigma
   !! build. Nothing here calls `mol%eris_packed` or `build_fock_direct`
   !! again: nothing in a Hessian-vector product needs a fresh AO integral
   !! pass once `a_block`/`b_block` and the density-independent inactive Fock
   !! exist, which is the whole point of building this state once rather than
   !! per application.
   !!
   !! **The four Hessian blocks, and where each one is cheap.**
   !!
   !! - Orbital-orbital: `one_index_fock` at the SA densities, unchanged --
   !!   exactly what `orbital_hessian` calls once per column, called here once
   !!   per input vector instead of once per unit vector.
   !! - CI-CI, block `J`: `2 w_J (H - E_J) x_J` through the existing
   !!   `sigma_vector`/`absorb_one_electron` machinery, no coupling between
   !!   states before projection.
   !! - Orbital-CI (into the orbital output, from `x_J`): the mixed partial
   !!   `d^2 E_SA / dkappa dx_J` equals `w_J` times the orbital-gradient
   !!   extraction of a generalised Fock built from the *symmetrised
   !!   transition density* between the reference `c_J` and `x_J`
   !!   (`transition_rdms`). Built here by `cheap_generalized_fock`,
   !!   which reuses `a_block`/`b_block`/the inactive Fock rather than
   !!   calling `generalized_fock` (which would redo the AO integral pass
   !!   every call) -- called with `delta_only = .true.`, because
   !!   `generalized_fock`'s inactive row is `2*(FI + FA(dm1))` and `FI` is a
   !!   density-independent constant that must NOT reappear when `(dm1, dm2)`
   !!   is a perturbation rather than a state (see that flag's docstring:
   !!   this was a real bug, caught by gate A comparing against
   !!   `orbital_hessian` on a converged SA point).
   !! - CI-orbital (into the CI output, from `kappa`): `2 w_J (H[kappa] -
   !!   <c_J|H[kappa]|c_J>) c_J`, with `H[kappa]` the active-space Hamiltonian
   !!   one-index-transformed along `kappa` -- `one_index_active_hamiltonian`
   !!   below, built from the same `a_block`/`b_block` slice `one_index_fock`
   !!   already differentiates, restricted to an active outer index instead of
   !!   a general one, plus the inactive Fock's own one-index transform (the
   !!   same construction `one_index_fock` uses for its `f_inactive`,
   !!   reproduced here rather than exposed as a byproduct of that routine).
   !!   Folded with the existing `absorb_one_electron` (linear, so folding the
   !!   *derivative* of `(h_eff, eri_act)` is the derivative of the folded
   !!   tensor) and applied with the existing `sigma_vector`.
   !!
   !! **Redundancy.** Every CI component, input and output, is projected onto
   !! the complement of `span{c_1, ..., c_N}` (`project_ci_block`, PySCF's
   !! `project_Aop`) and symmetrised under transpose (the singlet-restriction
   !! `davidson_lowest`'s `symmetrize_singlet` already imposes on every `c_J`;
   !! `c(ia,ib) -> (c(ia,ib)+c(ib,ia))/2` is that symmetrisation written
   !! directly on the `(na, nb)` array rather than through a flat vector).
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_casci, only: active_space_integrals
   use mqc_czt_mcscf, only: mcscf_fock_t, generalized_fock, orbital_gradient, &
                            rotation_parameters, rotation_matrix, &
                            one_index_fock, transformed_potential, sa_density_matrices
   use mqc_czt_mp2, only: transform_block
   use mqc_determinants, only: link_table_t, build_link_table
   use mqc_ci, only: absorb_one_electron, sigma_vector, ci_diagonal
   use mqc_rdm, only: active_space_rdms, transition_rdms
   implicit none
   private

   public :: sa_hessian_t
   public :: build_sa_hessian
   public :: destroy_sa_hessian
   public :: sa_hessian_n_param
   public :: sa_gradient
   public :: sa_hessian_apply
   public :: sa_hessian_precondition
   public :: project_ci_block
   public :: cheap_generalized_fock   !! Exposed for the tests, against `generalized_fock`
   public :: one_index_active_hamiltonian   !! Exposed for the tests, against `active_space_integrals`
   public :: gather_from_general
      !! Exposed for `mqc_czt_sa_nac`: the same orbital-gradient extraction it
      !! applies to `cheap_generalized_fock`'s output when building the
      !! interstate-coupling Z-vector's orbital right-hand side.

   real(dp), parameter :: CURVATURE_FLOOR = 1.0e-3_dp
      !! Smallest magnitude the diagonal preconditioner divides by, matching
      !! `mqc_orbital_rotation`'s `MIN_CURVATURE` in spirit: a mode softer than
      !! this is not trusted to precondition itself.

   type :: sa_hessian_t
      !! Everything a Hessian-vector product needs, built once at a converged
      !! SA-CASSCF point
      integer :: n_ao = 0, n_mo = 0
      integer :: n_inactive = 0, n_active = 0, n_alpha = 0, n_beta = 0
      integer :: n_states = 0, n_rot = 0, n_det = 0
      type(czt_molecule_t) :: mol
      real(dp), allocatable :: orbitals(:, :)         !! (n_ao, n_mo), the reference
      real(dp), allocatable :: weights(:)             !! (n_states)
      real(dp), allocatable :: ci_vectors(:, :, :)    !! (na, nb, n_states), normalised
      real(dp), allocatable :: active_energies(:)     !! (n_states), total minus core
      real(dp) :: core_energy = 0.0_dp
      type(link_table_t) :: alpha, beta
      integer, allocatable :: rows(:), cols(:)        !! From `rotation_parameters`
      type(mcscf_fock_t) :: fock_sa                   !! At the SA densities
      real(dp), allocatable :: dm1_sa(:, :), dm2_sa(:, :, :, :)
      real(dp), allocatable :: a_block(:, :, :, :), b_block(:, :, :, :)
         !! `(n_mo, n_occ, n_mo, n_occ)`, as `orbital_hessian` builds them
      real(dp), allocatable :: h_eff(:, :), eri_act(:, :, :, :)   !! The reference active Hamiltonian
      real(dp), allocatable :: folded(:, :)           !! From `absorb_one_electron`
      real(dp), allocatable :: diagonal(:, :)         !! From `ci_diagonal`, (na, nb)
      real(dp), allocatable :: orbital_hessian_diag(:)   !! (n_rot), approximate (see `build_sa_hessian`)
   end type sa_hessian_t

contains

   pure function sa_hessian_n_param(state) result(n)
      !! The flat parameter length, `n_rot + n_states * n_det`
      type(sa_hessian_t), intent(in) :: state
      integer :: n

      n = state%n_rot + state%n_states*state%n_det
   end function sa_hessian_n_param

   subroutine build_sa_hessian(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                               ci_vectors, energies, weights, state, error)
      !! Build every reusable piece of the SA Hessian-vector product
      !!
      !! `ci_vectors`/`energies` are a converged `run_czt_casscf` SA result's
      !! `ci_vectors`/`energies` (total, core plus active); `orbitals` its
      !! `orbitals`. Restricted to a complete active space and `n_alpha ==
      !! n_beta`, as `run_czt_casscf`'s own state-averaged path is.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)     !! (na, nb, >= n_states)
      real(dp), intent(in) :: energies(:)             !! (>= n_states), total
      real(dp), intent(in) :: weights(:)              !! (n_states)
      type(sa_hessian_t), intent(out) :: state
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: eri_packed(:, :)
      integer :: n_occ, j, l, p, q

      if (error%has_error()) return
      if (n_alpha /= n_beta) then
         call error%set(ERROR_VALIDATION, "sa_hessian: state averaging needs equal "// &
                        "active alpha and beta electrons ("//to_char(n_alpha)//" and "// &
                        to_char(n_beta)//" given).")
         return
      end if
      if (size(ci_vectors, 3) < size(weights) .or. size(energies) < size(weights)) then
         call error%set(ERROR_VALIDATION, "sa_hessian: fewer CI vectors or energies "// &
                        "than the "//to_char(size(weights))//" states requested.")
         return
      end if

      state%n_ao = size(orbitals, 1)
      state%n_mo = size(orbitals, 2)
      state%n_inactive = n_inactive
      state%n_active = n_active
      state%n_alpha = n_alpha
      state%n_beta = n_beta
      state%n_states = size(weights)
      n_occ = n_inactive + n_active

      state%mol = mol
      state%orbitals = orbitals
      state%weights = weights
      state%ci_vectors = ci_vectors(:, :, 1:state%n_states)

      call build_link_table(n_active, n_alpha, state%alpha, error)
      call build_link_table(n_active, n_beta, state%beta, error)
      if (error%has_error()) return
      state%n_det = state%alpha%n_strings*state%beta%n_strings

      call rotation_parameters(state%n_mo, n_inactive, n_active, state%rows, state%cols)
      state%n_rot = size(state%rows)

      call sa_density_matrices(state%ci_vectors, weights, state%alpha, state%beta, &
                               state%dm1_sa, state%dm2_sa, error)
      if (error%has_error()) return

      call generalized_fock(mol, orbitals, n_inactive, n_active, state%dm1_sa, &
                            state%dm2_sa, state%fock_sa, error)
      if (error%has_error()) return

      call mol%eris_packed(eri_packed)
      call transform_block(eri_packed, orbitals, orbitals(:, 1:n_occ), orbitals, &
                           orbitals(:, 1:n_occ), state%a_block)
      call transform_block(eri_packed, orbitals, orbitals, orbitals(:, 1:n_occ), &
                           orbitals(:, 1:n_occ), state%b_block)
      deallocate (eri_packed)

      call active_space_integrals(mol, orbitals, n_inactive, n_active, state%h_eff, &
                                  state%eri_act, state%core_energy, error)
      if (error%has_error()) return
      call absorb_one_electron(state%h_eff, state%eri_act, n_alpha + n_beta, &
                               state%folded, error)
      if (error%has_error()) return
      call ci_diagonal(state%h_eff, state%eri_act, state%alpha, state%beta, &
                       state%diagonal, error)
      if (error%has_error()) return

      allocate (state%active_energies(state%n_states))
      do j = 1, state%n_states
         state%active_energies(j) = energies(j) - state%core_energy
      end do

      ! The orbital preconditioner: the one-electron approximation to the
      ! Hessian diagonal, 2 (n_q f_pp + n_p f_qq) - 2 (F_pp + F_qq), with
      ! f = FI + FA and F the generalised Fock. For an inactive-virtual pair
      ! it is 4 (f_aa - f_ii). Building the exact diagonal would take one
      ! `one_index_fock` per rotation.
      allocate (state%orbital_hessian_diag(state%n_rot))
      do l = 1, state%n_rot
         p = state%rows(l)
         q = state%cols(l)
         state%orbital_hessian_diag(l) = 2.0_dp*(state%fock_sa%occupation(q)* &
                                                 (state%fock_sa%inactive(p, p) + state%fock_sa%active(p, p)) &
                                                 + state%fock_sa%occupation(p)* &
                                                 (state%fock_sa%inactive(q, q) + state%fock_sa%active(q, q))) &
                                         - 2.0_dp*(state%fock_sa%general(p, p) + state%fock_sa%general(q, q))
      end do
   end subroutine build_sa_hessian

   subroutine destroy_sa_hessian(state)
      !! Release the molecule and the excitation tables the state owns
      type(sa_hessian_t), intent(inout) :: state

      call state%mol%destroy()
      call state%alpha%destroy()
      call state%beta%destroy()
   end subroutine destroy_sa_hessian

   subroutine cheap_generalized_fock(state, dm1, dm2, general, delta_only)
      !! `fock%general` for an arbitrary `(dm1, dm2)` pair on the active
      !! space, reusing `a_block`/`b_block` and the density-independent
      !! inactive Fock instead of a fresh AO integral pass
      !!
      !! Exactly `generalized_fock`'s row assembly, with the active mean
      !! field read out of `transformed_potential(a_block, b_block, ...)`
      !! (already the MO-transformed integrals) rather than built by
      !! `build_fock_direct` on a fresh AO-basis density, and the inactive
      !! Fock taken from `state%fock_sa%inactive`, which never depends on the
      !! active density at all.
      !!
      !! **`generalized_fock` is affine in `(dm1, dm2)`, not linear**: the
      !! inactive row is `2*(FI + FA(dm1))`, and `FI` does not depend on the
      !! density it is fed at all, so it survives even at `dm1 = dm2 = 0`.
      !! That constant is exactly right for a genuine state's density (the
      !! default here, and what the consistency gate below checks against
      !! `generalized_fock` itself), but wrong for a *delta* density -- a
      !! transition density between a reference state and a CI trial vector,
      !! which is a first-order perturbation, not a state, and whose own
      !! derivative through the orbital-independent `FI` constant is exactly
      !! zero. `delta_only = .true.` drops that constant from the inactive
      !! row for exactly that case (the orbital-CI coupling block in
      !! `sa_hessian_apply_one`); the active row needs no such correction,
      !! since both its terms are already genuinely linear in `(dm1, dm2)`
      !! (a fixed matrix multiplying the density, never multiplying itself).
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: dm1(:, :)              !! (n_active, n_active)
      real(dp), intent(in) :: dm2(:, :, :, :)        !! (n_active^4), chemist order
      real(dp), allocatable, intent(out) :: general(:, :)   !! (n_mo, n_mo)
      logical, intent(in), optional :: delta_only
         !! `dm1`/`dm2` are a perturbation, not a state: omit the
         !! density-independent `2*FI` inactive-row constant. Default false.

      real(dp), allocatable :: occupations(:, :), potential(:, :)
      real(dp) :: accumulated
      logical :: skip_constant
      integer :: n_mo, n_occ, n_inactive, n_active, i, t, u, v, w, n

      skip_constant = .false.
      if (present(delta_only)) skip_constant = delta_only

      n_mo = state%n_mo
      n_inactive = state%n_inactive
      n_active = state%n_active
      n_occ = n_inactive + n_active

      allocate (occupations(n_mo, n_mo))
      occupations = 0.0_dp
      occupations(n_inactive + 1:n_occ, n_inactive + 1:n_occ) = dm1
      allocate (potential(n_mo, n_occ))
      call transformed_potential(state%a_block, state%b_block, n_occ, occupations, potential)

      allocate (general(n_mo, n_mo))
      general = 0.0_dp
      do i = 1, n_inactive
         do n = 1, n_mo
            if (skip_constant) then
               general(i, n) = 2.0_dp*potential(n, i)
            else
               general(i, n) = 2.0_dp*(state%fock_sa%inactive(n, i) + potential(n, i))
            end if
         end do
      end do

      do t = 1, n_active
         do n = 1, n_mo
            accumulated = 0.0_dp
            do u = 1, n_active
               accumulated = accumulated + dm1(t, u)*state%fock_sa%inactive(n, n_inactive + u)
            end do
            do w = 1, n_active
               do v = 1, n_active
                  do u = 1, n_active
                     accumulated = accumulated + dm2(t, u, v, w)* &
                                   state%a_block(n, n_inactive + u, n_inactive + v, n_inactive + w)
                  end do
               end do
            end do
            general(n_inactive + t, n) = accumulated
         end do
      end do

      deallocate (occupations, potential)
   end subroutine cheap_generalized_fock

   subroutine one_index_active_hamiltonian(state, kappa, dh_eff, deri_act)
      !! The active-space Hamiltonian's one- and two-electron integrals,
      !! differentiated along one orbital rotation `kappa`
      !!
      !! `dh_eff = f_inactive(kappa)_active`: `h_eff` (`active_space_integrals`)
      !! is already `C_active^T FI_AO C_active` -- the *whole* inactive Fock
      !! `h_ao + J(D_inactive) - K(D_inactive)/2`, not the bare integral with
      !! the mean field added separately -- so its one-index transform is
      !! exactly `one_index_fock`'s (unreturned) `f_inactive` intermediate,
      !! reproduced here: the commutator of `fock_sa%inactive` with `kappa`
      !! (which already carries `h_ao`'s own re-expression; `h_ao` does not
      !! depend on the inactive density, so nothing about it needs a second,
      !! separate commutator) plus `transformed_potential` of the rotated
      !! inactive density. A first version *did* differentiate a separate
      !! bare-`h_ao` term on top of this and double-counted it -- caught by
      !! `test_one_index_active_hamiltonian`'s finite-difference gate, which
      !! is why that check exists. `deri_act` is the same four-term
      !! construction `one_index_fock` uses for its `eri_gaaa`, with the
      !! general index restricted to the active range instead of every MO:
      !! `a_block`/`b_block` already hold `(n u|v w)` for `n` general, so
      !! slicing `n` to active costs nothing extra in integrals, only in how
      !! much of the same loop is kept.
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: kappa(:, :)            !! (n_mo, n_mo), antisymmetric
      real(dp), allocatable, intent(out) :: dh_eff(:, :)          !! (n_active, n_active)
      real(dp), allocatable, intent(out) :: deri_act(:, :, :, :)  !! (n_active^4)

      real(dp), allocatable :: d_inactive(:, :), v_inactive(:, :)
      real(dp), allocatable :: left(:, :), right(:, :), comm_fi(:, :)
      integer :: n_mo, n_occ, n_inactive, n_active, t, u, v, w, ta, ua, va, wa

      n_mo = state%n_mo
      n_inactive = state%n_inactive
      n_active = state%n_active
      n_occ = n_inactive + n_active

      ! The inactive density's own one-index transform, verbatim from
      ! `one_index_fock`: non-zero only where exactly one index is inactive.
      allocate (d_inactive(n_mo, n_mo))
      d_inactive = 0.0_dp
      d_inactive(:, 1:n_inactive) = 2.0_dp*kappa(:, 1:n_inactive)
      d_inactive(1:n_inactive, :) = d_inactive(1:n_inactive, :) - 2.0_dp*kappa(1:n_inactive, :)

      allocate (v_inactive(n_mo, n_occ))
      call transformed_potential(state%a_block, state%b_block, n_occ, d_inactive, v_inactive)

      allocate (left(n_mo, n_mo), right(n_mo, n_mo))
      call pic_gemm(state%fock_sa%inactive, kappa, left)
      call pic_gemm(kappa, state%fock_sa%inactive, right)
      allocate (comm_fi(n_mo, n_occ))
      comm_fi = left(:, 1:n_occ) - right(:, 1:n_occ) + v_inactive

      allocate (dh_eff(n_active, n_active))
      do u = 1, n_active
         ua = n_inactive + u
         do t = 1, n_active
            ta = n_inactive + t
            dh_eff(t, u) = comm_fi(ta, ua)
         end do
      end do

      allocate (deri_act(n_active, n_active, n_active, n_active))
      do w = 1, n_active
         wa = n_inactive + w
         do v = 1, n_active
            va = n_inactive + v
            do u = 1, n_active
               ua = n_inactive + u
               do t = 1, n_active
                  ta = n_inactive + t
                  deri_act(t, u, v, w) = &
                     dot_product(kappa(:, ta), state%a_block(:, ua, va, wa)) &
                     + dot_product(kappa(:, ua), state%b_block(ta, :, va, wa)) &
                     + dot_product(kappa(:, va), state%a_block(ta, ua, :, wa)) &
                     + dot_product(kappa(:, wa), state%a_block(ta, ua, :, va))
               end do
            end do
         end do
      end do

      deallocate (d_inactive, v_inactive, left, right, comm_fi)
   end subroutine one_index_active_hamiltonian

   subroutine project_ci_block(state, v)
      !! Remove the `span{c_1, ..., c_N}` component of `v`, in place
      !!
      !! PySCF's `project_Aop`: every reference state, not only the one `v`
      !! is attached to, since rotations among the averaged states are the
      !! redundant directions this removes. One pass, since the reference
      !! vectors are Davidson eigenvectors of a symmetric operator and
      !! orthonormal to working precision.
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(inout) :: v(:, :)     !! (na, nb)

      real(dp) :: overlap
      integer :: k

      do k = 1, state%n_states
         overlap = sum(state%ci_vectors(:, :, k)*v)
         v = v - overlap*state%ci_vectors(:, :, k)
      end do
   end subroutine project_ci_block

   subroutine unpack_kappa(state, flat, kappa)
      !! The antisymmetric `(n_mo, n_mo)` matrix a flat rotation vector packs
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: flat(:)        !! (n_rot)
      real(dp), allocatable, intent(out) :: kappa(:, :)

      integer :: l

      allocate (kappa(state%n_mo, state%n_mo))
      kappa = 0.0_dp
      do l = 1, state%n_rot
         kappa(state%rows(l), state%cols(l)) = flat(l)
         kappa(state%cols(l), state%rows(l)) = -flat(l)
      end do
   end subroutine unpack_kappa

   pure function gather_from_general(state, general) result(flat)
      !! `flat(l) = 2 (general(cols(l), rows(l)) - general(rows(l), cols(l)))`
      !!
      !! The same extraction `orbital_gradient` and `orbital_hessian`'s inner
      !! loop use, applied directly to a raw `(n_mo, n_mo)` matrix rather than
      !! through the `mcscf_fock_t` wrapper, so it can be reused on
      !! `cheap_generalized_fock`'s output as well as on `one_index_fock`'s.
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: general(:, :)
      real(dp), allocatable :: flat(:)

      integer :: l

      allocate (flat(state%n_rot))
      do l = 1, state%n_rot
         flat(l) = 2.0_dp*(general(state%cols(l), state%rows(l)) &
                           - general(state%rows(l), state%cols(l)))
      end do
   end function gather_from_general

   subroutine sa_gradient(state, kappa, x, gradient_orb, gradient_ci, error)
      !! `dE_SA/d(kappa, {x_J})` at an arbitrary point, not only a converged
      !! one -- what the finite-difference gate needs, since it evaluates
      !! this on both sides of the reference point
      !!
      !! `kappa` moves the orbitals as `C -> C exp(kappa)`; each `c_J + x_J`
      !! is renormalised before its density is built, so the reported energy
      !! and gradient are the genuine `<c'|H|c'>/<c'|c'>` at every point, not
      !! only to first order in `x_J` -- required for the Hessian (a second
      !! derivative) to come out right from finite differences of this first
      !! derivative. At `kappa = 0`, every `x_J = 0`, this must reproduce the
      !! reference point's (near-)vanishing SA gradient.
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: kappa(:, :)            !! (n_mo, n_mo), antisymmetric
      real(dp), intent(in) :: x(:, :, :)             !! (na, nb, n_states)
      real(dp), allocatable, intent(out) :: gradient_orb(:, :)   !! (n_mo, n_mo)
      real(dp), allocatable, intent(out) :: gradient_ci(:, :, :)  !! (na, nb, n_states)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: rotation(:, :), new_orbitals(:, :)
      real(dp), allocatable :: h_eff(:, :), eri_act(:, :, :, :), folded(:, :)
      real(dp), allocatable :: c_j(:, :), sigma_j(:, :), dm1_j(:, :), dm2_j(:, :, :, :)
      real(dp), allocatable :: dm1_sa(:, :), dm2_sa(:, :, :, :)
      type(mcscf_fock_t) :: fock
      real(dp) :: norm2, e_active, core_energy
      integer :: j, na, nb

      if (error%has_error()) return
      na = state%alpha%n_strings
      nb = state%beta%n_strings

      call rotation_matrix(kappa, rotation)
      allocate (new_orbitals(state%n_ao, state%n_mo))
      call pic_gemm(state%orbitals, rotation, new_orbitals)

      call active_space_integrals(state%mol, new_orbitals, state%n_inactive, state%n_active, &
                                  h_eff, eri_act, core_energy, error)
      if (error%has_error()) return
      call absorb_one_electron(h_eff, eri_act, state%n_alpha + state%n_beta, folded, error)
      if (error%has_error()) return

      allocate (gradient_ci(na, nb, state%n_states))
      allocate (sigma_j(na, nb))
      do j = 1, state%n_states
         c_j = state%ci_vectors(:, :, j) + x(:, :, j)
         norm2 = sum(c_j*c_j)
         call sigma_vector(folded, c_j, state%alpha, state%beta, sigma_j, error)
         if (error%has_error()) return
         e_active = sum(c_j*sigma_j)/norm2
         gradient_ci(:, :, j) = (2.0_dp*state%weights(j)/norm2)*(sigma_j - e_active*c_j)

         call active_space_rdms(c_j/sqrt(norm2), state%alpha, state%beta, dm1_j, dm2_j, error)
         if (error%has_error()) return
         if (j == 1) then
            dm1_sa = state%weights(1)*dm1_j
            dm2_sa = state%weights(1)*dm2_j
         else
            dm1_sa = dm1_sa + state%weights(j)*dm1_j
            dm2_sa = dm2_sa + state%weights(j)*dm2_j
         end if
      end do
      deallocate (sigma_j)

      call generalized_fock(state%mol, new_orbitals, state%n_inactive, state%n_active, &
                            dm1_sa, dm2_sa, fock, error)
      if (error%has_error()) return
      call orbital_gradient(fock, state%n_inactive, state%n_active, gradient_orb)

      deallocate (new_orbitals, rotation, h_eff, eri_act, folded, dm1_sa, dm2_sa)
   end subroutine sa_gradient

   subroutine sa_hessian_apply(state, x, hx, error)
      !! `H_SA` applied to a block of flat parameter vectors, one column at a
      !! time
      !!
      !! **The redundancy projection (`project_ci_block`, applied to every
      !! CI input and output here) is exact only for equal weights.** With
      !! unequal `w_J`, an in-space rotation between states `J` and `K` has
      !! curvature `(w_J - w_K)(E_K - E_J)` (`developer_sa_casscf.rst`,
      !! "Redundancy projection") and is not a true null direction of
      !! `E_SA`; projecting it out here still runs, but the resulting
      !! operator is then the Hessian of `E_SA` restricted to the subspace
      !! orthogonal to every reference state, not the unconstrained Hessian.
      !! The SA-CASSCF gradient is equal-weights only, matching PySCF's own
      !! refusal for that case; this routine itself does not
      !! refuse unequal weights, since nothing about the block formulas
      !! stops working, only the projection's interpretation changes.
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: x(:, :)     !! (n_param, n_vec)
      real(dp), intent(out) :: hx(:, :)   !! (n_param, n_vec)
      type(error_t), intent(inout) :: error

      integer :: iv

      do iv = 1, size(x, 2)
         call sa_hessian_apply_one(state, x(:, iv), hx(:, iv), error)
         if (error%has_error()) return
      end do
   end subroutine sa_hessian_apply

   subroutine sa_hessian_apply_one(state, xflat, hxflat, error)
      !! `H_SA` applied to one flat parameter vector
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: xflat(:)      !! (n_param)
      real(dp), intent(out) :: hxflat(:)    !! (n_param)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: kappa(:, :), transformed(:, :)
      real(dp), allocatable :: hx_orb(:), delta_orb(:)
      real(dp), allocatable :: x_j(:, :), hx_j(:, :)
      real(dp), allocatable :: tdm1a(:, :), tdm2a(:, :, :, :)
      real(dp), allocatable :: tdm1(:, :), tdm2(:, :, :, :)
      real(dp), allocatable :: dh_eff(:, :), deri_act(:, :, :, :), folded_delta(:, :)
      real(dp), allocatable :: sigma_j(:, :), sigma_delta(:, :)
      real(dp), allocatable :: fock_transition(:, :)
      real(dp) :: e_delta
      integer :: na, nb, j, seg0, seg1

      if (error%has_error()) return
      na = state%alpha%n_strings
      nb = state%beta%n_strings

      call unpack_kappa(state, xflat(1:state%n_rot), kappa)

      ! ---- orbital-orbital: one_index_fock at the SA densities -----------
      ! `transformed` is a plain assumed-shape `intent(out)` in `one_index_fock`
      ! (not allocatable there), so it must already have the right shape.
      allocate (transformed(state%n_mo, state%n_mo))
      call one_index_fock(state%a_block, state%b_block, state%fock_sa, state%dm1_sa, &
                          state%dm2_sa, state%n_inactive, state%n_active, kappa, transformed)
      hx_orb = gather_from_general(state, transformed)

      ! ---- the kappa-differentiated active Hamiltonian, shared by every
      ! state's CI-orbital block -------------------------------------------
      call one_index_active_hamiltonian(state, kappa, dh_eff, deri_act)
      call absorb_one_electron(dh_eff, deri_act, state%n_alpha + state%n_beta, &
                               folded_delta, error)
      if (error%has_error()) return

      allocate (hx_j(na, nb), sigma_j(na, nb), sigma_delta(na, nb))
      do j = 1, state%n_states
         seg0 = state%n_rot + (j - 1)*state%n_det + 1
         seg1 = state%n_rot + j*state%n_det
         x_j = reshape(xflat(seg0:seg1), [na, nb])
         call project_ci_block(state, x_j)
         x_j = 0.5_dp*(x_j + transpose(x_j))

         ! ---- CI-CI: 2 w_J (H - E_J) x_J, block-diagonal across states ---
         call sigma_vector(state%folded, x_j, state%alpha, state%beta, sigma_j, error)
         if (error%has_error()) return
         hx_j = 2.0_dp*state%weights(j)*(sigma_j - state%active_energies(j)*x_j)

         ! ---- CI-orbital: 2 w_J (H[kappa] - <c_J|H[kappa]|c_J>) c_J ------
         call sigma_vector(folded_delta, state%ci_vectors(:, :, j), state%alpha, state%beta, &
                           sigma_delta, error)
         if (error%has_error()) return
         e_delta = sum(state%ci_vectors(:, :, j)*sigma_delta)
         hx_j = hx_j + 2.0_dp*state%weights(j)*(sigma_delta - e_delta*state%ci_vectors(:, :, j))

         call project_ci_block(state, hx_j)
         hx_j = 0.5_dp*(hx_j + transpose(hx_j))
         hxflat(seg0:seg1) = reshape(hx_j, [state%n_det])

         ! ---- orbital-CI: w_J * gradient-extraction of the generalised
         ! Fock built from the symmetrised c_J/x_J transition density ------
         ! The (x_J, c_J) ordering is the (c_J, x_J) one transposed:
         ! dm1 -> dm1^T, dm2(p,q,r,s) -> dm2(q,p,s,r).
         call transition_rdms(state%ci_vectors(:, :, j), x_j, state%alpha, state%beta, &
                              tdm1a, tdm2a, error)
         if (error%has_error()) return
         tdm1 = tdm1a + transpose(tdm1a)
         tdm2 = tdm2a + reshape(tdm2a, shape(tdm2a), order=[2, 1, 4, 3])
         call cheap_generalized_fock(state, tdm1, tdm2, fock_transition, delta_only=.true.)
         delta_orb = gather_from_general(state, fock_transition)
         hx_orb = hx_orb + state%weights(j)*delta_orb
      end do

      hxflat(1:state%n_rot) = hx_orb

      deallocate (kappa, transformed, hx_orb, dh_eff, deri_act, folded_delta)
      deallocate (hx_j, sigma_j, sigma_delta, x_j)
      deallocate (tdm1a, tdm2a, tdm1, tdm2, fock_transition, delta_orb)
   end subroutine sa_hessian_apply_one

   subroutine sa_hessian_precondition(state, x, px)
      !! A diagonal preconditioner: the orbital block from the one-electron
      !! approximation to the orbital Hessian's diagonal, the
      !! CI block from `2 w_J (H_diag - E_J)`, both floored so a near-zero
      !! curvature is not divided by
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: x(:, :)      !! (n_param, n_vec)
      real(dp), intent(out) :: px(:, :)    !! (n_param, n_vec)

      real(dp) :: d
      integer :: l, j, iv, seg0, det, ia, ib, na

      na = state%alpha%n_strings
      do iv = 1, size(x, 2)
         do l = 1, state%n_rot
            d = state%orbital_hessian_diag(l)
            if (abs(d) < CURVATURE_FLOOR) d = sign(CURVATURE_FLOOR, d)
            px(l, iv) = x(l, iv)/d
         end do
         do j = 1, state%n_states
            seg0 = state%n_rot + (j - 1)*state%n_det
            do det = 1, state%n_det
               ia = mod(det - 1, na) + 1
               ib = (det - 1)/na + 1
               d = 2.0_dp*state%weights(j)*(state%diagonal(ia, ib) - state%active_energies(j))
               if (abs(d) < CURVATURE_FLOOR) d = sign(CURVATURE_FLOOR, d)
               px(seg0 + det, iv) = x(seg0 + det, iv)/d
            end do
         end do
      end do
   end subroutine sa_hessian_precondition

end module mqc_czt_sa_hessian
