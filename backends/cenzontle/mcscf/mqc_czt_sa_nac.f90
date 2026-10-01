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
   !! convention -- see `sa_casscf_nac_pair`'s comment on the point). It is
   !! what `use_etfs`
   !! (electron translation factors) omits; PySCF includes it by default, so
   !! `include_csf` defaults true here to match. Unlike `h_base` and the
   !! response terms, it is not translationally invariant on its own (an
   !! overlap-derivative contraction between two different orbital sets, not
   !! a Hellmann-Feynman-type energy derivative) -- expected, not a bug; see
   !! `developer_sa_casscf.rst`.
   !!
   !! **What is fused and what is not.** `czt_sa_casscf_nacs` builds the SA
   !! Hessian state once (`build_sa_hessian`) and shares it across every
   !! requested pair, the same expensive precompute
   !! `czt_sa_casscf_gradients` shares across roots. Each pair still solves
   !! its own single-column Z-vector equation rather than joining a block
   !! solve, and its own two derivative-integral sweeps (`h_base` and the two
   !! response terms) rather than a cross-pair fused sweep: the gradient's
   !! fusion (`czt_sa_casscf_gradients`) batches structurally identical right-hand sides across
   !! roots; a NAC's right-hand side and its `h_base` both depend on the pair
   !! `(I,J)` itself in a way that would need re-deriving that fusion from
   !! scratch. Left for a later pass if NACs for many pairs turn out to
   !! dominate a run's cost.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t, atom_ao_blocks
   use mqc_czt_gradient, only: one_electron_deriv, DERIV_OVLP
   use mqc_czt_mcscf_gradient, only: active_two_electron_gradient
   use mqc_czt_sa_hessian, only: sa_hessian_t, build_sa_hessian, destroy_sa_hessian, &
                                 sa_hessian_n_param, cheap_generalized_fock, &
                                 project_ci_block, gather_from_general
   use mqc_czt_sa_gradient, only: sa_zvector_solve, response_separable_gradient, &
                                  orbital_response_gradient, ci_response_gradient, &
                                  UNEQUAL_WEIGHT_TOL
   use mqc_ci, only: sigma_vector
   use mqc_rdm, only: transition_rdms
   implicit none
   private

   ! `czt_sa_casscf_nac` is public for the tests, which check one pair at a
   ! time against PySCF; the bridge calls `czt_sa_casscf_nacs`.
   public :: czt_sa_casscf_nac
   public :: czt_sa_casscf_nacs

   real(dp), parameter :: DEFAULT_CG_TOL = 1.0e-10_dp
      !! Relative residual each pair's Z-vector solve stops at, as
      !! `mqc_czt_sa_gradient`'s
   integer, parameter :: DEFAULT_CG_MAX_ITER = 200
      !! As `mqc_czt_sa_gradient`'s
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
      !! states of a converged SA-CASSCF, building its own SA Hessian state --
      !! `czt_sa_casscf_nacs` shares one across several pairs instead.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)          !! (n_ao, n_mo), the SA orbitals
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      real(dp), intent(in) :: ci_vectors(:, :, :)     !! (na, nb, >= n_states)
      real(dp), intent(in) :: energies(:)             !! (>= n_states), total
      real(dp), intent(in) :: weights(:)              !! (n_states)
      integer, intent(in) :: state_i, state_j         !! 1-based, `state_i /= state_j`
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

      type(sa_hessian_t) :: state

      if (error%has_error()) return
      call build_sa_hessian(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                            ci_vectors, energies, weights, state, error)
      if (error%has_error()) return

      call sa_casscf_nac_pair(state, orbitals, n_inactive, n_active, weights, state_i, &
                              state_j, energies, coupling, interstate_coupling, csf_term, &
                              energy_difference, error, include_csf, cg_tol, cg_max_iter, &
                              cg_iterations, cg_residual)
      call destroy_sa_hessian(state)
   end subroutine czt_sa_casscf_nac

   subroutine czt_sa_casscf_nacs(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                 ci_vectors, energies, weights, pairs, couplings, &
                                 interstate_couplings, csf_terms, energy_differences, error, &
                                 include_csf, cg_tol, cg_max_iter, cg_iterations, cg_residual)
      !! Every pair in `pairs`, sharing one SA Hessian state -- see the module
      !! docstring for what is and is not fused across pairs
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

      type(sa_hessian_t) :: state
      real(dp), allocatable :: coupling(:, :), interstate(:, :), csf(:, :)
      integer :: n_states, n_pairs, ip, its
      real(dp) :: res, ediff

      if (error%has_error()) return
      n_states = size(weights)
      n_pairs = size(pairs, 2)
      if (present(cg_iterations)) cg_iterations = 0
      if (present(cg_residual)) cg_residual = 0.0_dp

      do ip = 1, n_pairs
         if (any(pairs(:, ip) < 1) .or. any(pairs(:, ip) > n_states)) then
            call error%set(ERROR_VALIDATION, "sa_casscf_nacs: pair ("// &
                           to_char(pairs(1, ip))//","//to_char(pairs(2, ip))// &
                           ") is not within the "//to_char(n_states)//" averaged states.")
            return
         end if
      end do

      call build_sa_hessian(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                            ci_vectors, energies, weights, state, error)
      if (error%has_error()) return

      allocate (couplings(3, mol%natm, n_pairs))
      allocate (interstate_couplings(3, mol%natm, n_pairs))
      allocate (csf_terms(3, mol%natm, n_pairs))
      allocate (energy_differences(n_pairs))

      do ip = 1, n_pairs
         call sa_casscf_nac_pair(state, orbitals, n_inactive, n_active, weights, &
                                 pairs(1, ip), pairs(2, ip), energies, coupling, interstate, &
                                 csf, ediff, error, include_csf, cg_tol, cg_max_iter, its, res)
         if (error%has_error()) then
            call destroy_sa_hessian(state)
            return
         end if
         couplings(:, :, ip) = coupling
         interstate_couplings(:, :, ip) = interstate
         csf_terms(:, :, ip) = csf
         energy_differences(ip) = ediff
         if (present(cg_iterations)) cg_iterations(ip) = its
         if (present(cg_residual)) cg_residual(ip) = res
      end do

      call destroy_sa_hessian(state)
   end subroutine czt_sa_casscf_nacs

   subroutine sa_casscf_nac_pair(state, orbitals, n_inactive, n_active, weights, state_i, &
                                 state_j, energies, coupling, interstate_coupling, csf_term, &
                                 energy_difference, error, include_csf, cg_tol, cg_max_iter, &
                                 cg_iterations, cg_residual)
      !! The Z-vector machinery on an already-built `state`, shared by
      !! `czt_sa_casscf_nac` (which builds one just for this call) and
      !! `czt_sa_casscf_nacs` (which shares one across every pair)
      type(sa_hessian_t), intent(in) :: state
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: state_i, state_j
      real(dp), intent(in) :: energies(:)
      real(dp), allocatable, intent(out) :: coupling(:, :)
      real(dp), allocatable, intent(out) :: interstate_coupling(:, :)
      real(dp), allocatable, intent(out) :: csf_term(:, :)
      real(dp), intent(out) :: energy_difference
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: include_csf
      real(dp), intent(in), optional :: cg_tol
      integer, intent(in), optional :: cg_max_iter
      integer, intent(out), optional :: cg_iterations
      real(dp), intent(out), optional :: cg_residual

      real(dp), allocatable :: tdm1a(:, :), tdm2a(:, :, :, :)
      real(dp), allocatable :: tdm1_sym(:, :), tdm2_sym(:, :, :, :)
      real(dp), allocatable :: c_active(:, :), c_core(:, :)
      real(dp), allocatable :: d_core(:, :), d_active_sa(:, :)
      real(dp), allocatable :: fock_transition(:, :), grad_flat(:)
      real(dp), allocatable :: rhs(:, :), x(:, :)
      real(dp), allocatable :: kappa_bar(:, :), xbar(:, :, :)
      real(dp), allocatable :: sigma_i(:, :), sigma_j(:, :)
      real(dp), allocatable :: orb_response(:, :), ci_response(:, :)
      real(dp), allocatable :: work(:, :)
      real(dp), allocatable :: tm1_antisym(:, :), tm1_ao(:, :)
      real(dp), allocatable :: s1(:, :, :)
      integer, allocatable :: offsets(:), counts(:)
      real(dp) :: use_tol, residual
      logical :: want_csf
      integer :: n_states, use_max_iter, iterations
      integer :: n_ao, n_mo, na, nb, n_param, n_occ, l, j, seg0, seg1
      integer :: iatom, comp, p0, p1

      if (error%has_error()) return
      n_states = size(weights)
      if (present(cg_iterations)) cg_iterations = 0
      if (present(cg_residual)) cg_residual = 0.0_dp
      energy_difference = 0.0_dp

      if (state_i < 1 .or. state_i > n_states .or. state_j < 1 .or. state_j > n_states) then
         call error%set(ERROR_VALIDATION, "sa_casscf_nac: states "//to_char(state_i)// &
                        " and "//to_char(state_j)//" must both be among the "// &
                        to_char(n_states)//" averaged states.")
         return
      end if
      if (maxval(weights) - minval(weights) > UNEQUAL_WEIGHT_TOL) then
         call error%set(ERROR_VALIDATION, "sa_casscf_nac: unequal SA weights -- the "// &
                        "redundancy projection in the SA Hessian is exact only for "// &
                        "equal weights, so a nonadiabatic coupling is refused rather "// &
                        "than built from the wrong Lagrangian.")
         return
      end if

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_occ = n_inactive + n_active
      na = state%alpha%n_strings
      nb = state%beta%n_strings

      allocate (coupling(3, state%mol%natm), interstate_coupling(3, state%mol%natm))
      allocate (csf_term(3, state%mol%natm))
      coupling = 0.0_dp
      interstate_coupling = 0.0_dp
      csf_term = 0.0_dp
      energy_difference = energies(state_j) - energies(state_i)

      if (state_i == state_j) return   ! d_II is identically zero; PySCF's kernel() short-circuit

      want_csf = .true.
      if (present(include_csf)) want_csf = include_csf
      ! Refused before any work: `coupling` divides by this. At a conical
      ! intersection `h_IJ` is still finite, but returning it alone would make
      ! `coupling` mean something different pair by pair.
      if (abs(energy_difference) < DEGENERATE_ENERGY_TOL) then
         call error%set(ERROR_VALIDATION, "nonadiabatic coupling: states "// &
                        to_char(state_i)//" and "//to_char(state_j)//" are degenerate "// &
                        "(|E_J - E_I| below "//to_char(DEGENERATE_ENERGY_TOL)// &
                        " Hartree), where d_IJ = h_IJ/(E_J - E_I) is not defined.")
         return
      end if

      use_tol = DEFAULT_CG_TOL
      if (present(cg_tol)) use_tol = cg_tol
      use_max_iter = DEFAULT_CG_MAX_ITER
      if (present(cg_max_iter)) use_max_iter = cg_max_iter

      ! ---- the symmetrised transition density between I and J, PySCF's
      ! `make_fcasscf_nacs` convention: 1/2 (D_IJ + D_JI) ---------------------
      call transition_rdms(state%ci_vectors(:, :, state_i), state%ci_vectors(:, :, state_j), &
                           state%alpha, state%beta, tdm1a, tdm2a, error)
      if (error%has_error()) return
      tdm1_sym = 0.5_dp*(tdm1a + transpose(tdm1a))
      tdm2_sym = 0.5_dp*(tdm2a + reshape(tdm2a, shape(tdm2a), order=[2, 1, 4, 3]))

      allocate (c_active(n_ao, n_active))
      c_active = orbitals(:, n_inactive + 1:n_occ)
      allocate (d_core(n_ao, n_ao))
      d_core = 0.0_dp
      if (n_inactive > 0) then
         allocate (c_core(n_ao, n_inactive))
         c_core = orbitals(:, 1:n_inactive)
         call pic_gemm(c_core, c_core, d_core, transb="T", alpha=2.0_dp, beta=0.0_dp)
      end if
      allocate (d_active_sa(n_ao, n_ao), work(n_ao, n_active))
      call pic_gemm(c_active, state%dm1_sa, work, beta=0.0_dp)
      call pic_gemm(work, c_active, d_active_sa, transb="T", beta=0.0_dp)
      deallocate (work)

      ! ---- h_base: <I|dH/dR|J>, exactly the CI-response machinery's
      ! separable and active-active pieces, fed the I/J transition density --
      call response_separable_gradient(state%mol, d_core, d_active_sa, &
                                       build_d_ci_active(c_active, tdm1_sym), &
                                       build_weighted_transition(state, orbitals, tdm1_sym, &
                                                                 tdm2_sym), &
                                       interstate_coupling, error)
      if (error%has_error()) return
      call active_two_electron_gradient(state%mol, c_active, tdm2_sym, interstate_coupling, &
                                        error)
      if (error%has_error()) return

      ! ---- the Z-vector right-hand side: orbital block from the same
      ! transition-density generalised Fock, CI blocks from H_active acting on
      ! the OTHER state, swapped between I's and J's slots (see module
      ! docstring) ------------------------------------------------------------
      call cheap_generalized_fock(state, tdm1_sym, tdm2_sym, fock_transition, &
                                  delta_only=.true.)
      grad_flat = gather_from_general(state, fock_transition)

      n_param = sa_hessian_n_param(state)
      allocate (rhs(n_param, 1))
      rhs = 0.0_dp
      do l = 1, state%n_rot
         rhs(l, 1) = -grad_flat(l)
      end do

      allocate (sigma_j(na, nb), sigma_i(na, nb))
      call sigma_vector(state%folded, state%ci_vectors(:, :, state_j), state%alpha, &
                        state%beta, sigma_j, error)
      if (error%has_error()) return
      call project_ci_block(state, sigma_j)
      seg0 = state%n_rot + (state_i - 1)*state%n_det + 1
      seg1 = state%n_rot + state_i*state%n_det
      rhs(seg0:seg1, 1) = -reshape(sigma_j, [state%n_det])

      call sigma_vector(state%folded, state%ci_vectors(:, :, state_i), state%alpha, &
                        state%beta, sigma_i, error)
      if (error%has_error()) return
      call project_ci_block(state, sigma_i)
      seg0 = state%n_rot + (state_j - 1)*state%n_det + 1
      seg1 = state%n_rot + state_j*state%n_det
      rhs(seg0:seg1, 1) = -reshape(sigma_i, [state%n_det])

      call sa_zvector_solve(state, rhs, x, iterations, residual, use_tol, use_max_iter, error)
      if (present(cg_iterations)) cg_iterations = iterations
      if (present(cg_residual)) cg_residual = residual
      if (error%has_error()) return

      allocate (kappa_bar(n_mo, n_mo))
      kappa_bar = 0.0_dp
      do l = 1, state%n_rot
         kappa_bar(state%rows(l), state%cols(l)) = x(l, 1)
         kappa_bar(state%cols(l), state%rows(l)) = -x(l, 1)
      end do
      allocate (xbar(na, nb, n_states))
      do j = 1, n_states
         seg0 = state%n_rot + (j - 1)*state%n_det + 1
         seg1 = state%n_rot + j*state%n_det
         xbar(:, :, j) = reshape(x(seg0:seg1, 1), [na, nb])
      end do

      call orbital_response_gradient(state%mol, orbitals, n_inactive, n_active, state, &
                                     kappa_bar, orb_response, error)
      if (error%has_error()) return
      call ci_response_gradient(state%mol, orbitals, n_inactive, n_active, state, xbar, &
                                weights, ci_response, error)
      if (error%has_error()) return

      interstate_coupling = interstate_coupling + orb_response + ci_response

      ! ---- the CSF term: the antisymmetric part of the raw (unsymmetrised)
      ! transition 1-RDM against the AO overlap derivative, scaled by the
      ! energy gap, as PySCF's `nac_csf`. The antisymmetric part is taken as
      ! `transpose - plain` because this module's h_IJ = (E_J - E_I) d_IJ has
      ! the opposite overall sign to PySCF's raw h_IJ, which scales by
      ! E_I - E_J; the other pieces share the Z-vector machinery and carry
      ! the sign already. ----------------------------------------------------
      allocate (tm1_antisym(n_active, n_active))
      tm1_antisym = transpose(tdm1a) - tdm1a
      allocate (tm1_ao(n_ao, n_ao), work(n_ao, n_active))
      call pic_gemm(c_active, tm1_antisym, work, beta=0.0_dp)
      call pic_gemm(work, c_active, tm1_ao, transb="T", beta=0.0_dp)
      deallocate (work)

      call one_electron_deriv(state%mol, s1, DERIV_OVLP)
      s1 = -s1   ! the same sign `response_separable_assemble`'s Pulay term uses
      allocate (offsets(state%mol%natm), counts(state%mol%natm))
      call atom_ao_blocks(state%mol, offsets, counts)
      do iatom = 1, state%mol%natm
         p0 = offsets(iatom) + 1
         p1 = offsets(iatom) + counts(iatom)
         if (counts(iatom) == 0) cycle
         do comp = 1, 3
            csf_term(comp, iatom) = 0.5_dp*energy_difference* &
                                    sum(s1(p0:p1, :, comp)*tm1_ao(p0:p1, :))
         end do
      end do
      deallocate (offsets, counts, s1)

      if (want_csf) interstate_coupling = interstate_coupling + csf_term

      coupling = interstate_coupling/energy_difference
   end subroutine sa_casscf_nac_pair

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
