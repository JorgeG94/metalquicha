!! Nonadiabatic (derivative) couplings between two SA-CASSCF states
module mqc_czt_sa_nac
   !!
   !! `d_IJ = <Psi_I | d/dR Psi_J>` for two states `I /= J` of a converged
   !! SA-CASSCF, following the Lengsfield/Yarkony Lagrangian (Snyder,
   !! Hohenstein, Luehr & Martinez, J. Chem. Phys. 143, 154107 (2015); PySCF's
   !! `pyscf.nac.sacasscf`) and matching that reference's conventions exactly,
   !! confirmed by comparison against it:
   !!
   !!     (E_J - E_I) d_IJ = h_IJ = <I| dH/dR |J> + orbital response
   !!                              + CI response + CSF term
   !!
   !! **The interstate coupling `h_IJ`.** The bra-ket matrix element `<I|H|J>`
   !! restricted to the active space (the core/inactive part drops out: `I`
   !! and `J` are orthogonal determinant expansions on the same inactive core,
   !! so the core-Hamiltonian constant that multiplies `<I|J>` in a genuine
   !! state's energy multiplies zero here) is exactly the quantity
   !! `mqc_czt_sa_gradient`'s CI-response machinery already differentiates,
   !! with the caller-supplied "transition density between a Lagrange
   !! multiplier and a reference state" replaced by the transition density
   !! **between the two states themselves**, symmetrised the same way PySCF's
   !! `make_fcasscf_nacs` does: `dm1_sym = 1/2 (D_IJ + D_JI)`, `dm2_sym`
   !! likewise, with `D_IJ(p,q) = <I|E_pq|J>` (`transition_rdms`). Because
   !! `<I|H|J>` is already numerically the same whichever of `(D_IJ, D_JI)`
   !! computes it (a real, symmetric one/two-electron Hamiltonian sandwiched
   !! between two real vectors gives `<I|H|J> = <J|H|I>`), this symmetrisation
   !! changes nothing about the contracted energy, only makes the matrix
   !! `cheap_generalized_fock`/`gather_from_general` need (built for a
   !! Hermitian density) well-formed. `h_base` is then exactly
   !! `response_separable_gradient` (core-active and active-active-separable
   !! Coulomb/Pulay pieces, no `d_core_resp`: a transition density lives on the
   !! active space alone) plus `active_two_electron_gradient` on `dm2_sym`
   !! directly (not a cumulant: a transition density is not a state's own
   !! two-particle density).
   !!
   !! **Orbital and CI response.** The same Z-vector equation
   !! `H_SA [kappa_bar; xbar_1..xbar_N] = -RHS` the gradient uses
   !! (`sa_zvector_solve`, on `mqc_czt_sa_hessian`'s shared operator), but with
   !! a different right-hand side, because a different quantity is being
   !! differentiated:
   !!
   !!   - orbital block: `-gather_from_general(cheap_generalized_fock(dm1_sym,
   !!     dm2_sym, delta_only=.true.))`, the same extraction
   !!     `orbital_response_pieces` applies to a state's own generalised Fock,
   !!     here to the transition-density one (`delta_only` because, as with
   !!     the CI-response term, `<I|H|J>` has none of a genuine state's
   !!     density-independent inactive-Fock constant: it multiplies `<I|J> =
   !!     0`).
   !!   - CI block, state `I`'s slot: `-(H_active c_J)`, projected off every
   !!     averaged state (`project_ci_block`) -- `<I|H|J>` is *linear* in
   !!     `c_I`, so its gradient wrt `c_I` is simply `H_active c_J` (no factor
   !!     of 2: that factor belongs to differentiating a *quadratic* form
   !!     `<c|H|c>` wrt its one shared argument, which is not what this is).
   !!   - CI block, state `J`'s slot: `-(H_active c_I)`, by the same argument
   !!     with `I`/`J` swapped (`<I|H|J> = <J|H|I>` for a real Hamiltonian).
   !!   - Every other state's CI block: zero (`<I|H|J>` does not depend on it).
   !!
   !! Given `(kappa_bar, {xbar_J})`, the response pieces themselves are
   !! `mqc_czt_sa_gradient`'s `orbital_response_gradient` and
   !! `ci_response_gradient` **unchanged**: the Lagrangian's constraint terms
   !! `kappa_bar . dE_SA/dkappa` and `sum_J w_J xbar_J . (H-E_J)c_J` enforce
   !! the same SA stationarity condition regardless of what quantity the
   !! multipliers were solved for, so nothing about differentiating them
   !! depends on whether the right-hand side came from a root's own gradient
   !! or from this module's transition density.
   !!
   !! **The CSF term.** PySCF's `nac_csf`: the antisymmetric part of the *raw*
   !! (unsymmetrised) transition 1-RDM, `D_IJ^T - D_IJ`, dotted with the AO
   !! overlap derivative and scaled by `E_J - E_I` (that sign, not PySCF's own
   !! internal `E_I - E_J`, is what keeps this term consistent with `h_base`'s
   !! convention -- see `czt_sa_casscf_gradients_nacs`'s comment on the point). It is
   !! what `use_etfs`
   !! (electron translation factors) omits; PySCF includes it by default, so
   !! `include_csf` defaults true here to match. Unlike `h_base` and the
   !! response terms, it is not translationally invariant on its own (an
   !! overlap-derivative contraction between two different orbital sets, not
   !! a Hellmann-Feynman-type energy derivative) -- expected, not a bug; see
   !! `developer_sa_casscf.rst`.
   !!
   !! **Everything is fused.** `czt_sa_casscf_gradients_nacs` builds each
   !! pair's right-hand side and `<I|dH/dR|J>` densities without touching a
   !! derivative integral (`nac_pair_inputs`) and hands them to
   !! `sa_gradients_on_state` as extra Lagrangian columns beside the roots'
   !! gradients: one SA Hessian state, one block Z-vector solve, one
   !! separable derivative-integral sweep and one active two-body sweep for
   !! every root and every pair. Only the CSF term is added afterwards.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks
   use mqc_czt_gradient, only: one_electron_deriv, DERIV_OVLP
   use mqc_czt_sa_hessian, only: sa_hessian_t, build_sa_hessian, destroy_sa_hessian, &
                                 sa_hessian_n_param, cheap_generalized_fock, &
                                 project_ci_block, gather_from_general
   use mqc_czt_sa_gradient, only: sa_gradients_on_state, UNEQUAL_WEIGHT_TOL
   use mqc_ci, only: sigma_vector
   use mqc_rdm, only: transition_rdms
   implicit none
   private

   ! `czt_sa_casscf_nac` is public for the tests, which check one pair at a
   ! time against PySCF; the bridge calls `czt_sa_casscf_nacs`.
   public :: czt_sa_casscf_nac
   public :: czt_sa_casscf_nacs
   public :: czt_sa_casscf_gradients_nacs

   real(dp), parameter :: DEGENERATE_ENERGY_TOL = 1.0e-8_dp
      !! Hartree. `|E_J - E_I|` below this is refused: `d_IJ = h_IJ/(E_J - E_I)`
      !! would be Inf or NaN, or numerically meaningless, at a degeneracy

contains

   subroutine czt_sa_casscf_nac(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                ci_vectors, energies, weights, state_i, state_j, &
                                coupling, interstate_coupling, csf_term, energy_difference, &
                                error, include_csf, cg_tol, cg_max_iter, cg_iterations, &
                                cg_residual)
      !! The full NAC `d_IJ` (`coupling`) and its pieces between one pair of
      !! states of a converged SA-CASSCF: `czt_sa_casscf_nacs` on one pair
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)          !! (n_ao, n_mo), the SA orbitals
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)     !! (na, nb, >= n_states)
      real(dp), intent(in) :: energies(:)             !! (>= n_states), total
      real(dp), intent(in) :: weights(:)              !! (n_states)
      integer, intent(in) :: state_i, state_j         !! 1-based
      real(dp), allocatable, intent(out) :: coupling(:, :)             !! (3, n_atoms), d_IJ, 1/Bohr
      real(dp), allocatable, intent(out) :: interstate_coupling(:, :)  !! (3, n_atoms), h_IJ, Hartree/Bohr
      real(dp), allocatable, intent(out) :: csf_term(:, :)             !! (3, n_atoms), already energy-scaled
      real(dp), intent(out) :: energy_difference      !! E_J - E_I, Hartree
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: include_csf     !! Default true, matching PySCF's `use_etfs=False`
      real(dp), intent(in), optional :: cg_tol
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations
      real(dp), intent(out), optional :: cg_residual

      real(dp), allocatable :: couplings(:, :, :), interstates(:, :, :), csfs(:, :, :)
      real(dp), allocatable :: ediffs(:)
      integer :: its(1)
      real(dp) :: res(1)

      energy_difference = 0.0_dp
      if (present(cg_iterations)) cg_iterations = 0
      if (present(cg_residual)) cg_residual = 0.0_dp
      if (error%has_error()) return
      call czt_sa_casscf_nacs(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                              ci_vectors, energies, weights, reshape([state_i, state_j], [2, 1]), &
                              couplings, interstates, csfs, ediffs, error, include_csf, &
                              cg_tol, cg_max_iter, its, res)
      if (error%has_error()) return
      coupling = couplings(:, :, 1)
      interstate_coupling = interstates(:, :, 1)
      csf_term = csfs(:, :, 1)
      energy_difference = ediffs(1)
      if (present(cg_iterations)) cg_iterations = its(1)
      if (present(cg_residual)) cg_residual = res(1)
   end subroutine czt_sa_casscf_nac

   subroutine czt_sa_casscf_nacs(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                 ci_vectors, energies, weights, pairs, couplings, &
                                 interstate_couplings, csf_terms, energy_differences, error, &
                                 include_csf, cg_tol, cg_max_iter, cg_iterations, cg_residual)
      !! Every pair in `pairs` and no root gradients: `czt_sa_casscf_gradients_nacs`
      !! with no roots
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)
      real(dp), intent(in) :: energies(:)
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: pairs(:, :)               !! (2, n_pairs), 1-based [state_i, state_j]
      real(dp), allocatable, intent(out) :: couplings(:, :, :)             !! (3, natm, n_pairs)
      real(dp), allocatable, intent(out) :: interstate_couplings(:, :, :)  !! (3, natm, n_pairs)
      real(dp), allocatable, intent(out) :: csf_terms(:, :, :)             !! (3, natm, n_pairs)
      real(dp), allocatable, intent(out) :: energy_differences(:)          !! (n_pairs)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: include_csf
      real(dp), intent(in), optional :: cg_tol
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations(:)   !! (n_pairs)
      real(dp), intent(out), optional :: cg_residual(:)     !! (n_pairs)

      real(dp), allocatable :: gradients(:, :, :)
      integer :: no_roots(0)

      call czt_sa_casscf_gradients_nacs(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                        ci_vectors, energies, weights, no_roots, pairs, &
                                        gradients, couplings, interstate_couplings, csf_terms, &
                                        energy_differences, error, include_csf, cg_tol, &
                                        cg_max_iter, nac_cg_iterations=cg_iterations, &
                                        nac_cg_residual=cg_residual)
   end subroutine czt_sa_casscf_nacs

   subroutine czt_sa_casscf_gradients_nacs(mol, orbitals, n_inactive, n_active, n_alpha, &
                                           n_beta, ci_vectors, energies, weights, roots, &
                                           pairs, gradients, couplings, interstate_couplings, &
                                           csf_terms, energy_differences, error, include_csf, &
                                           cg_tol, cg_max_iter, cg_iterations, cg_residual, &
                                           nac_cg_iterations, nac_cg_residual)
      !! The gradients of every root in `roots` and the couplings of every
      !! pair in `pairs`, relaxed together by `sa_gradients_on_state`: one SA
      !! Hessian state, one block Z-vector solve, one separable
      !! derivative-integral sweep and one active two-body sweep for all of
      !! them. A pair with `state_i == state_j` is zero and takes no column.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)
      real(dp), intent(in) :: energies(:)
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: roots(:)                   !! 1-based
      integer, intent(in) :: pairs(:, :)                !! (2, n_pairs), 1-based
      real(dp), allocatable, intent(out) :: gradients(:, :, :)             !! (3, natm, size(roots))
      real(dp), allocatable, intent(out) :: couplings(:, :, :)             !! (3, natm, n_pairs)
      real(dp), allocatable, intent(out) :: interstate_couplings(:, :, :)  !! (3, natm, n_pairs)
      real(dp), allocatable, intent(out) :: csf_terms(:, :, :)             !! (3, natm, n_pairs)
      real(dp), allocatable, intent(out) :: energy_differences(:)          !! (n_pairs), E_J - E_I
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: include_csf
      real(dp), intent(in), optional :: cg_tol
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations(:)       !! (size(roots))
      real(dp), intent(out), optional :: cg_residual(:)         !! (size(roots))
      integer, intent(out), optional :: nac_cg_iterations(:)   !! (n_pairs)
      real(dp), intent(out), optional :: nac_cg_residual(:)     !! (n_pairs)

      type(sa_hessian_t) :: state
      real(dp), allocatable :: extra_rhs(:, :), extra_gamma(:, :, :, :, :)
      real(dp), allocatable :: extra_d(:, :, :), extra_w(:, :, :), tm1_ao(:, :, :)
      real(dp), allocatable :: extra_out(:, :, :), s1(:, :, :)
      integer, allocatable :: column(:), its(:), extra_states(:, :)
      real(dp), allocatable :: res(:)
      integer, allocatable :: offsets(:), counts(:)
      logical :: want_csf
      integer :: n_states, n_pairs, n_roots, n_extra, n_ao, n_param, ip, ie, natm
      integer :: iatom, comp, p0, p1

      if (error%has_error()) return
      n_states = size(weights)
      n_pairs = size(pairs, 2)
      n_roots = size(roots)
      natm = mol%natm
      n_ao = size(orbitals, 1)
      want_csf = .true.
      if (present(include_csf)) want_csf = include_csf
      if (present(cg_iterations)) cg_iterations = 0
      if (present(cg_residual)) cg_residual = 0.0_dp
      if (present(nac_cg_iterations)) nac_cg_iterations = 0
      if (present(nac_cg_residual)) nac_cg_residual = 0.0_dp

      do ip = 1, n_pairs
         if (any(pairs(:, ip) < 1) .or. any(pairs(:, ip) > n_states)) then
            call error%set(ERROR_VALIDATION, "sa_casscf_nacs: pair ("// &
                           to_char(pairs(1, ip))//","//to_char(pairs(2, ip))// &
                           ") is not within the "//to_char(n_states)//" averaged states.")
            return
         end if
      end do
      if (maxval(weights) - minval(weights) > UNEQUAL_WEIGHT_TOL) then
         call error%set(ERROR_VALIDATION, "sa_casscf_nac: unequal SA weights -- the "// &
                        "redundancy projection in the SA Hessian is exact only for "// &
                        "equal weights, so a nonadiabatic coupling is refused rather "// &
                        "than built from the wrong Lagrangian.")
         return
      end if

      allocate (couplings(3, natm, n_pairs), interstate_couplings(3, natm, n_pairs))
      allocate (csf_terms(3, natm, n_pairs), energy_differences(n_pairs))
      couplings = 0.0_dp
      interstate_couplings = 0.0_dp
      csf_terms = 0.0_dp
      do ip = 1, n_pairs
         energy_differences(ip) = energies(pairs(2, ip)) - energies(pairs(1, ip))
      end do

      ! Refused before any work: `couplings` divides by these. At a conical
      ! intersection `h_IJ` is still finite, but returning it alone would make
      ! `couplings` mean something different pair by pair.
      do ip = 1, n_pairs
         if (pairs(1, ip) == pairs(2, ip)) cycle
         if (abs(energy_differences(ip)) < DEGENERATE_ENERGY_TOL) then
            call error%set(ERROR_VALIDATION, "nonadiabatic coupling: states "// &
                           to_char(pairs(1, ip))//" and "//to_char(pairs(2, ip))// &
                           " are degenerate (|E_J - E_I| below "// &
                           to_char(DEGENERATE_ENERGY_TOL)//" Hartree), where d_IJ = "// &
                           "h_IJ/(E_J - E_I) is not defined.")
            return
         end if
      end do

      ! d_II is identically zero (PySCF's kernel() short-circuit): no column.
      allocate (column(n_pairs))
      n_extra = 0
      do ip = 1, n_pairs
         column(ip) = 0
         if (pairs(1, ip) /= pairs(2, ip)) then
            n_extra = n_extra + 1
            column(ip) = n_extra
         end if
      end do

      call build_sa_hessian(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                            ci_vectors, energies, weights, state, error)
      if (error%has_error()) return

      n_param = sa_hessian_n_param(state)
      allocate (extra_rhs(n_param, n_extra))
      allocate (extra_gamma(n_active, n_active, n_active, n_active, n_extra))
      allocate (extra_d(n_ao, n_ao, n_extra), extra_w(n_ao, n_ao, n_extra))
      allocate (tm1_ao(n_ao, n_ao, n_extra), extra_states(2, n_extra))
      do ip = 1, n_pairs
         ie = column(ip)
         if (ie == 0) cycle
         extra_states(:, ie) = pairs(:, ip)
         call nac_pair_inputs(state, orbitals, n_inactive, n_active, pairs(1, ip), &
                              pairs(2, ip), extra_rhs(:, ie), extra_gamma(:, :, :, :, ie), &
                              extra_d(:, :, ie), extra_w(:, :, ie), tm1_ao(:, :, ie), error)
         if (error%has_error()) then
            call destroy_sa_hessian(state)
            return
         end if
      end do

      allocate (its(n_roots + n_extra), res(n_roots + n_extra))
      call sa_gradients_on_state(state, orbitals, n_inactive, n_active, weights, roots, &
                                 gradients, error, cg_tol, cg_max_iter, its, res, &
                                 extra_rhs=extra_rhs, extra_gamma=extra_gamma, &
                                 extra_d_active=extra_d, extra_weighted=extra_w, &
                                 extra_out=extra_out, extra_states=extra_states)
      call destroy_sa_hessian(state)
      if (present(cg_iterations)) cg_iterations = its(1:n_roots)
      if (present(cg_residual)) cg_residual = res(1:n_roots)
      if (error%has_error()) return

      ! ---- the CSF term: the antisymmetric part of the raw (unsymmetrised)
      ! transition 1-RDM against the AO overlap derivative, scaled by the
      ! energy gap, as PySCF's `nac_csf`. The antisymmetric part is taken as
      ! `transpose - plain` because this module's h_IJ = (E_J - E_I) d_IJ has
      ! the opposite overall sign to PySCF's raw h_IJ, which scales by
      ! E_I - E_J; the other pieces share the Z-vector machinery and carry
      ! the sign already. ----------------------------------------------------
      if (n_extra > 0) then
         call one_electron_deriv(mol, s1, DERIV_OVLP)
         s1 = -s1   ! the same sign `response_separable_assemble`'s Pulay term uses
         allocate (offsets(natm), counts(natm))
         call atom_ao_blocks(mol, offsets, counts)
      end if
      do ip = 1, n_pairs
         ie = column(ip)
         if (ie == 0) cycle
         if (present(nac_cg_iterations)) nac_cg_iterations(ip) = its(n_roots + ie)
         if (present(nac_cg_residual)) nac_cg_residual(ip) = res(n_roots + ie)
         do iatom = 1, natm
            if (counts(iatom) == 0) cycle
            p0 = offsets(iatom) + 1
            p1 = offsets(iatom) + counts(iatom)
            do comp = 1, 3
               csf_terms(comp, iatom, ip) = 0.5_dp*energy_differences(ip)* &
                                            sum(s1(p0:p1, :, comp)*tm1_ao(p0:p1, :, ie))
            end do
         end do
         interstate_couplings(:, :, ip) = extra_out(:, :, ie)
         if (want_csf) interstate_couplings(:, :, ip) = interstate_couplings(:, :, ip) &
                                                        + csf_terms(:, :, ip)
         couplings(:, :, ip) = interstate_couplings(:, :, ip)/energy_differences(ip)
      end do
   end subroutine czt_sa_casscf_gradients_nacs

   subroutine nac_pair_inputs(state, orbitals, n_inactive, n_active, state_i, state_j, rhs, &
                              gamma_base, d_base, weighted_base, tm1_ao, error)
      !! Pair `(I, J)`'s Lagrangian column for `sa_gradients_on_state`, built
      !! without a derivative integral: the Z-vector right-hand side, the
      !! `<I|dH/dR|J>` base terms, and the AO antisymmetric transition density
      !! the CSF term needs. See the module docstring for each.
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, state_i, state_j
      real(dp), intent(out) :: rhs(:)                     !! (n_param)
      real(dp), intent(out) :: gamma_base(:, :, :, :)     !! (n_active^4), 1/2 (D_IJ + D_JI)
      real(dp), intent(out) :: d_base(:, :)               !! (n_ao, n_ao)
      real(dp), intent(out) :: weighted_base(:, :)        !! (n_ao, n_ao)
      real(dp), intent(out) :: tm1_ao(:, :)               !! (n_ao, n_ao), C (D_IJ^T - D_IJ) C^T
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: tdm1a(:, :), tdm2a(:, :, :, :), tdm1_sym(:, :)
      real(dp), allocatable :: c_active(:, :), work(:, :), grad_flat(:), fock_transition(:, :)
      real(dp), allocatable :: sigma(:, :), tm1_antisym(:, :)
      integer :: n_ao, na, nb, seg0, seg1

      n_ao = size(orbitals, 1)
      na = state%alpha%n_strings
      nb = state%beta%n_strings

      ! PySCF's `make_fcasscf_nacs` convention: 1/2 (D_IJ + D_JI).
      call transition_rdms(state%ci_vectors(:, :, state_i), state%ci_vectors(:, :, state_j), &
                           state%alpha, state%beta, tdm1a, tdm2a, error)
      if (error%has_error()) return
      tdm1_sym = 0.5_dp*(tdm1a + transpose(tdm1a))
      gamma_base = 0.5_dp*(tdm2a + reshape(tdm2a, shape(tdm2a), order=[2, 1, 4, 3]))

      c_active = orbitals(:, n_inactive + 1:n_inactive + n_active)
      d_base = build_d_ci_active(c_active, tdm1_sym)
      weighted_base = build_weighted_transition(state, orbitals, tdm1_sym, gamma_base)

      ! Orbital block from the transition-density generalised Fock; CI blocks
      ! from H_active acting on the OTHER state, swapped between I's and J's
      ! slots.
      call cheap_generalized_fock(state, tdm1_sym, gamma_base, fock_transition, &
                                  delta_only=.true.)
      grad_flat = gather_from_general(state, fock_transition)
      rhs = 0.0_dp
      rhs(1:state%n_rot) = -grad_flat(1:state%n_rot)

      allocate (sigma(na, nb))
      call sigma_vector(state%folded, state%ci_vectors(:, :, state_j), state%alpha, &
                        state%beta, sigma, error)
      if (error%has_error()) return
      call project_ci_block(state, sigma)
      seg0 = state%n_rot + (state_i - 1)*state%n_det + 1
      seg1 = state%n_rot + state_i*state%n_det
      rhs(seg0:seg1) = -reshape(sigma, [state%n_det])

      call sigma_vector(state%folded, state%ci_vectors(:, :, state_i), state%alpha, &
                        state%beta, sigma, error)
      if (error%has_error()) return
      call project_ci_block(state, sigma)
      seg0 = state%n_rot + (state_j - 1)*state%n_det + 1
      seg1 = state%n_rot + state_j*state%n_det
      rhs(seg0:seg1) = -reshape(sigma, [state%n_det])

      tm1_antisym = transpose(tdm1a) - tdm1a
      allocate (work(n_ao, n_active))
      call pic_gemm(c_active, tm1_antisym, work, beta=0.0_dp)
      call pic_gemm(work, c_active, tm1_ao, transb="T", beta=0.0_dp)
   end subroutine nac_pair_inputs

   function build_d_ci_active(c_active, tdm1_sym) result(d_ci_active)
      !! The AO transition density `C tdm1_sym C^T`, on the active block alone
      real(dp), intent(in) :: c_active(:, :)      !! (n_ao, n_active)
      real(dp), intent(in) :: tdm1_sym(:, :)       !! (n_active, n_active)
      real(dp), allocatable :: d_ci_active(:, :)   !! (n_ao, n_ao)

      real(dp), allocatable :: work(:, :)

      allocate (work(size(c_active, 1), size(c_active, 2)))
      allocate (d_ci_active(size(c_active, 1), size(c_active, 1)))
      call pic_gemm(c_active, tdm1_sym, work, beta=0.0_dp)
      call pic_gemm(work, c_active, d_ci_active, transb="T", beta=0.0_dp)
      deallocate (work)
   end function build_d_ci_active

   function build_weighted_transition(state, orbitals, tdm1_sym, tdm2_sym) result(weighted)
      !! The AO energy-weighted (overlap/Pulay) matrix for the `h_base` term:
      !! `C sym(cheap_generalized_fock(tdm1_sym, tdm2_sym)) C^T`, exactly
      !! `ci_response_pieces`'s `weighted_ci`, built from the I/J transition
      !! density instead of the multiplier/state one
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: orbitals(:, :)
      real(dp), intent(in) :: tdm1_sym(:, :)
      real(dp), intent(in) :: tdm2_sym(:, :, :, :)
      real(dp), allocatable :: weighted(:, :)   !! (n_ao, n_ao)

      real(dp), allocatable :: fock(:, :), sym_fock(:, :), work(:, :)
      integer :: n_ao, n_mo

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      call cheap_generalized_fock(state, tdm1_sym, tdm2_sym, fock, delta_only=.true.)
      allocate (sym_fock(n_mo, n_mo))
      sym_fock = 0.5_dp*(fock + transpose(fock))
      allocate (weighted(n_ao, n_ao), work(n_ao, n_mo))
      call pic_gemm(orbitals, sym_fock, work, beta=0.0_dp)
      call pic_gemm(work, orbitals, weighted, transb="T", beta=0.0_dp)
      deallocate (work)
   end function build_weighted_transition

end module mqc_czt_sa_nac
