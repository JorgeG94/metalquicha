!! The generalised Fock matrix, the orbital gradient and the orbital Hessian
module mqc_czt_mcscf
   !! What CASSCF adds to CASCI is that the orbitals move. CASCI minimises the
   !! energy over the CI coefficients at fixed orbitals; the derivative with
   !! respect to the orbitals is generally not zero, and driving it to zero is
   !! the whole of the extra work.
   !!
   !! Orbital changes are parametrised as `C -> C exp(kappa)` with `kappa`
   !! antisymmetric, so the transformation is orthogonal for any `kappa` and the
   !! orbitals cannot drift out of orthonormality however large a step is taken.
   !! The derivative of the energy with respect to `kappa_pq` is
   !!
   !!     g_pq = 2 (F_qp - F_pq)
   !!
   !! where `F` is the generalised Fock matrix. **That sign is fixed by finite
   !! differences, not asserted**, in `test_mqc_mcscf.f90`: both orderings
   !! appear in the literature because both `C exp(kappa)` and `C exp(-kappa)`
   !! are used as the parametrisation, and getting it backwards gives an
   !! optimiser that climbs. The rows of `F` are built differently depending on
   !! what the orbital is:
   !!
   !!     F_in = 2 (FI_ni + FA_ni)                            n inactive
   !!     F_tn = sum_u D_tu FI_nu + sum_uvw d_tuvw (nu|vw)    t active
   !!     F_an = 0                                            a virtual
   !!
   !! with `FI` the inactive Fock -- the closed-shell Fock of the doubly
   !! occupied orbitals -- and `FA` the mean field of the active density. A
   !! virtual orbital is empty in every determinant and contributes nothing, so
   !! its row vanishes; that is why the gradient has no virtual-virtual block
   !! and why rotations among virtual orbitals are redundant.
   !!
   !! **Most orbital rotations are redundant and must be left out.** Mixing two
   !! inactive orbitals does not change a single determinant, and neither does
   !! mixing two active orbitals in a *complete* active space, because the CI
   !! spans every distribution of the electrons over them either way. Including
   !! redundant parameters does not merely waste effort: the Hessian is singular
   !! along them, so a Newton step is undefined and even a gradient step wanders
   !! along directions the energy does not depend on. The three non-redundant
   !! blocks are inactive-active, inactive-virtual and active-virtual.
   !!
   !! **The Hessian is the same expression again, not a second one.** Expanding
   !! `C exp(kappa)` shows the term of the energy quadratic in `kappa` to be the
   !! term linear in it, evaluated on integrals differentiated once, so
   !! `one_index_fock` is `generalized_fock` with each orbital in turn replaced
   !! by `C kappa` and the Hessian is read out of it through the gradient
   !! formula above. Nothing here can have a Hessian that disagrees with its
   !! gradient.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_gemm, pic_dgemm_x
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_czt_integrals, only: czt_molecule_t
   use mqc_czt_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_czt_mp2, only: transform_block
   use mqc_czt_casci, only: run_czt_casci, casci_result_t, &
                            run_czt_ormas_ci
   use mqc_determinants, only: link_table_t, build_link_table
   use mqc_rdm, only: active_space_rdms, spin_squared
   use mqc_ormas_space, only: ormas_space_t, build_ormas_space
   use mqc_ormas_ci, only: ormas_density_matrices
   use pic_logger, only: logger => global_logger
   use pic_lapack_interfaces, only: pic_syev
   use mqc_orbital_rotation, only: rotation_matrix, level_shifted_step, &
                                   rotation_hessian_t, subspace_newton_step, &
                                   MIN_CURVATURE, SADDLE_CURVATURE, MAX_ROTATION, &
                                   MIN_ROTATION, ENERGY_RESOLUTION, TRUST_GROWTH
   implicit none
   private

   public :: mcscf_fock_t
   public :: generalized_fock
   public :: orbital_gradient
   public :: rotation_matrix
   public :: is_redundant
   public :: subspace_of
   public :: rotation_parameters
   public :: orbital_hessian
   public :: orbital_hessian_from_blocks
   public :: mo_integral_blocks
   public :: fock_from_blocks
   public :: orbital_hessian_operator_t
   public :: approximate_hessian_diagonal
   public :: iterative_newton_step
   public :: run_czt_casscf
   public :: casscf_result_t
   public :: natural_orbitals
   ! Public for `mqc_czt_sa_hessian`, which reuses rather than reimplements
   ! them: its orbital-orbital block is exactly `one_index_fock`, and its
   ! one-index-transformed active Hamiltonian needs the same inactive-density
   ! potential, from `transformed_potential`.
   public :: one_index_fock
   public :: transformed_potential
   public :: transformed_potential_many
   public :: one_index_fock_many
   public :: sa_density_matrices

   ! The step-control constants and the matrix exponential are
   ! `mqc_orbital_rotation`'s, used from there rather than declared here: the
   ! second-order SCF takes the same trust-region Newton step on the same
   ! parametrisation, and two copies of `MIN_CURVATURE` would be two things to
   ! keep equal. The numbers are unchanged; see that module for what each is
   ! for.

   type :: casscf_result_t
      !! What an orbital optimisation leaves behind
      real(dp) :: energy = 0.0_dp
      real(dp) :: core_energy = 0.0_dp
      real(dp) :: active_energy = 0.0_dp
      real(dp) :: gradient_norm = 0.0_dp     !! Largest element, at exit
      real(dp), allocatable :: orbitals(:, :)
      real(dp), allocatable :: ci_vector(:, :)
         !! (n_alpha_strings, n_beta_strings). Left unallocated by a restricted
         !! space, whose determinants are not a rectangle -- `ci_flat` carries
         !! it there instead, exactly as in `casci_result_t`.
      real(dp), allocatable :: ci_flat(:)
         !! (n_determinants), the vector of a restricted space
      real(dp), allocatable :: dm1(:, :)     !! Active one-particle density
      real(dp), allocatable :: dm2(:, :, :, :)
         !! Active two-particle density, spin-traced. Consistent with `dm1`,
         !! `orbitals` and `energy` when the optimisation converged, because the
         !! loop tests the gradient and leaves before touching the orbitals. On
         !! a run that ran out of iterations instead, both densities are one
         !! orbital step behind -- consistent with each other, which is what a
         !! cumulant needs, but not with the orbitals they are reported beside.
         !!
         !! Under state averaging (`n_states > 1`, below) these are the
         !! **state-averaged** densities `D_SA = sum_J w_J D_J`, `d_SA = sum_J
         !! w_J d_J` -- what the orbital optimiser actually used -- and not any
         !! one root's own density. A per-root density is not carried here;
         !! `energy` is likewise `E_SA`, not root 1's energy.
      real(dp), allocatable :: energies(:)
         !! Every state's total energy, ascending, when `n_states > 1`;
         !! `weights(J)` belongs to `energies(J)`. Unallocated otherwise.
      real(dp), allocatable :: ci_vectors(:, :, :)
         !! (n_alpha_strings, n_beta_strings, n_states): every state's CI
         !! vector at the final orbitals, when `n_states > 1`
      real(dp), allocatable :: spins(:)
         !! `<S^2>` of each state in `energies`, from its CI vector
         !! (`spin_squared`). Same shape as `energies`.
      integer :: iterations = 0
      integer :: n_determinants = 0
      logical :: converged = .false.
      logical :: stalled = .false.
         !! The optimisation stopped because no step downhill could be found,
         !! rather than because it ran out of iterations -- so more iterations
         !! will not help. A step worth less than `ENERGY_RESOLUTION` is taken
         !! untested, so this means the surface defeated the quadratic model
         !! and not that the gradient outran what the energy can resolve.
   end type casscf_result_t

   type :: mcscf_fock_t
      !! The generalised Fock matrix and the pieces it was built from
      real(dp), allocatable :: general(:, :)    !! (n_mo, n_mo), `F_mn`
      real(dp), allocatable :: inactive(:, :)   !! (n_mo, n_mo), `FI` in the MO basis
      real(dp), allocatable :: active(:, :)     !! (n_mo, n_mo), `FA` in the MO basis
      real(dp), allocatable :: occupation(:)
         !! (n_mo). Two for inactive, `D_tt` for active, zero for virtual. Worth
         !! having explicitly because it is what makes a rotation redundant when
         !! two orbitals share it.
   end type mcscf_fock_t

   type, extends(rotation_hessian_t) :: orbital_hessian_operator_t
      !! The orbital Hessian at fixed CI, applied to vectors without being built
      !!
      !! A vector holds one value per non-redundant rotation, in `rows`/`cols`
      !! order, and its image is that vector times `orbital_hessian_from_blocks`'s
      !! columns before that routine averages the matrix with its transpose.
      !! The integral blocks are the caller's, pointed to rather than copied.
      real(dp), pointer, contiguous :: a_block(:, :, :, :) => null()
         !! (n_mo, n_occ, n_mo, n_occ): `(p q|r s)`, `q` and `s` occupied
      real(dp), pointer, contiguous :: b_block(:, :, :, :) => null()
         !! (n_mo, n_mo, n_occ, n_occ): `(p q|r s)`, `r` and `s` occupied
      type(mcscf_fock_t) :: fock
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      integer :: n_inactive = 0
      integer :: n_active = 0
      integer, allocatable :: rows(:), cols(:)
   contains
      procedure :: apply => orbital_hessian_apply
   end type orbital_hessian_operator_t

   integer, parameter :: ITERATIVE_HESSIAN_ABOVE = 800
      !! Rotation count above which `run_czt_casscf` takes its Newton step
      !! iteratively instead of building and diagonalising the Hessian
   integer, parameter :: ITERATIVE_SUBSPACE = 40
      !! Largest Krylov subspace of an iterative Newton step
   real(dp), parameter :: ITERATIVE_TOLERANCE = 1.0e-4_dp
      !! Residual of the projected Newton equations, relative to the gradient,
      !! at which an iterative step stops expanding

contains

   pure function is_redundant(p, q, n_inactive, n_active, subspaces) result(redundant)
      !! Whether rotating `p` into `q` changes the wave function at all
      !!
      !! It does not if both are inactive or both virtual: the wave function is
      !! built from those sets, not from the orbitals inside them.
      !!
      !! **The active-active case depends on the space.** A complete active
      !! space distributes its electrons over its orbitals in every way there
      !! is, so mixing two of them reaches nothing new and the rotation is
      !! redundant. Restrict the occupations and that stops being true the
      !! moment the two orbitals fall in different subspaces -- the wave
      !! function then does distinguish them, and the rotation is a real
      !! parameter with a real gradient. Two orbitals *within* one subspace are
      !! redundant again, because a subspace is complete in itself.
      !!
      !! Getting this wrong is not loud: treat a real parameter as redundant and
      !! the optimiser stops somewhere that is not a stationary point, and treat
      !! a redundant one as real and the Hessian acquires a null direction.
      integer, intent(in) :: p, q, n_inactive, n_active
      integer, intent(in), optional :: subspaces(:)
         !! Active orbital each subspace starts at, ascending, as
         !! `keywords.mcscf.ormas.subspaces` gives it. Absent is one subspace
         !! covering everything.
      logical :: redundant

      integer :: class_p, class_q

      class_p = orbital_class(p, n_inactive, n_active)
      class_q = orbital_class(q, n_inactive, n_active)
      redundant = (class_p == class_q)

      if (redundant .and. class_p == 2 .and. present(subspaces)) then
         redundant = subspace_of(p - n_inactive, subspaces) == &
            subspace_of(q - n_inactive, subspaces)
      end if
   end function is_redundant

   pure function subspace_of(active_orbital, subspaces) result(which)
      !! Which subspace an active orbital belongs to, counting from 1
      !!
      !! `subspaces` is ascending and its first entry is 1, so the answer is the
      !! last entry not past the orbital.
      integer, intent(in) :: active_orbital
      integer, intent(in) :: subspaces(:)
      integer :: which

      integer :: k

      which = 1
      do k = 1, size(subspaces)
         if (subspaces(k) <= active_orbital) which = k
      end do
   end function subspace_of

   pure function orbital_class(p, n_inactive, n_active) result(class_index)
      !! 1 inactive, 2 active, 3 virtual
      integer, intent(in) :: p, n_inactive, n_active
      integer :: class_index

      if (p <= n_inactive) then
         class_index = 1
      else if (p <= n_inactive + n_active) then
         class_index = 2
      else
         class_index = 3
      end if
   end function orbital_class

   subroutine generalized_fock(mol, orbitals, n_inactive, n_active, dm1, dm2, &
                               fock, error)
      !! Build `F_mn`, and the inactive and active Fock matrices with it
      ! TODO(mqc): calls `mol%eris_packed` for its (n u|v w) block and makes two
      ! direct AO Fock builds. A caller that also needs the orbital Hessian
      ! should use `mo_integral_blocks` and `fock_from_blocks` instead, as
      ! `run_czt_casscf` does.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)      !! (n_ao, n_mo)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1(:, :)           !! Active one-particle density
      real(dp), intent(in) :: dm2(:, :, :, :)     !! Active two-particle density
      type(mcscf_fock_t), intent(out) :: fock
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: h_ao(:, :), bounds(:, :)
      real(dp), allocatable :: d_inactive(:, :), d_active(:, :)
      real(dp), allocatable :: f_inactive_ao(:, :), f_active_ao(:, :), zero_h(:, :)
      real(dp), allocatable :: c_inactive(:, :), c_active(:, :), work(:, :)
      real(dp), allocatable :: eri_gaaa(:, :, :, :), eri_packed(:, :)
      type(direct_stats_t) :: stats
      real(dp) :: accumulated
      integer :: n_ao, n_mo, i, t, u, v, w, n

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)

      if (n_inactive + n_active > n_mo) then
         call error%set(ERROR_VALIDATION, to_char(n_inactive)//" inactive plus "// &
                        to_char(n_active)//" active orbitals is more than the "// &
                        to_char(n_mo)//" available.")
         return
      end if

      allocate (c_active(n_ao, n_active))
      c_active = orbitals(:, n_inactive + 1:n_inactive + n_active)

      call mol%core_hamiltonian(h_ao)
      call schwarz_bounds(mol, bounds, error)
      if (error%has_error()) return

      ! The inactive Fock: the closed-shell Fock of the doubly occupied
      ! orbitals, which `build_fock_direct` computes as H + J - K/2.
      allocate (d_inactive(n_ao, n_ao), f_inactive_ao(n_ao, n_ao))
      d_inactive = 0.0_dp
      if (n_inactive > 0) then
         allocate (c_inactive(n_ao, n_inactive))
         c_inactive = orbitals(:, 1:n_inactive)
         call pic_gemm(c_inactive, c_inactive, d_inactive, transb="T", &
                       alpha=2.0_dp, beta=0.0_dp)
         deallocate (c_inactive)
      end if
      call build_fock_direct(mol, h_ao, d_inactive, bounds, f_inactive_ao, stats, error)
      if (error%has_error()) return

      ! The active mean field: the same J - K/2 built from the active density
      ! and with no one-electron part, which is what the zero core Hamiltonian
      ! leaves out.
      allocate (d_active(n_ao, n_ao), f_active_ao(n_ao, n_ao), zero_h(n_ao, n_ao))
      allocate (work(n_ao, n_active))
      call pic_gemm(c_active, dm1, work)
      call pic_gemm(work, c_active, d_active, transb="T")
      zero_h = 0.0_dp
      call build_fock_direct(mol, zero_h, d_active, bounds, f_active_ao, stats, error)
      if (error%has_error()) return
      deallocate (work)

      ! Both into the molecular orbital basis, over every orbital: the gradient
      ! couples occupied orbitals to virtual ones, so the virtual columns are
      ! needed even though virtual rows are zero.
      allocate (fock%inactive(n_mo, n_mo), fock%active(n_mo, n_mo))
      allocate (work(n_ao, n_mo))
      call pic_gemm(f_inactive_ao, orbitals, work)
      call pic_gemm(orbitals, work, fock%inactive, transa="T")
      call pic_gemm(f_active_ao, orbitals, work)
      call pic_gemm(orbitals, work, fock%active, transa="T")
      deallocate (work)

      ! (nu|vw): one general index, three active.
      call mol%eris_packed(eri_packed)
      call transform_block(eri_packed, orbitals, c_active, c_active, c_active, eri_gaaa)

      allocate (fock%general(n_mo, n_mo), fock%occupation(n_mo))
      fock%general = 0.0_dp
      fock%occupation = 0.0_dp

      do i = 1, n_inactive
         fock%occupation(i) = 2.0_dp
         do n = 1, n_mo
            fock%general(i, n) = 2.0_dp*(fock%inactive(n, i) + fock%active(n, i))
         end do
      end do

      do t = 1, n_active
         fock%occupation(n_inactive + t) = dm1(t, t)
         do n = 1, n_mo
            accumulated = 0.0_dp
            do u = 1, n_active
               accumulated = accumulated + dm1(t, u)*fock%inactive(n, n_inactive + u)
            end do
            do w = 1, n_active
               do v = 1, n_active
                  do u = 1, n_active
                     accumulated = accumulated + dm2(t, u, v, w)*eri_gaaa(n, u, v, w)
                  end do
               end do
            end do
            fock%general(n_inactive + t, n) = accumulated
         end do
      end do

      deallocate (h_ao, bounds, d_inactive, d_active, f_inactive_ao, f_active_ao)
      deallocate (zero_h, c_active, eri_gaaa, eri_packed)
   end subroutine generalized_fock

   subroutine orbital_gradient(fock, n_inactive, n_active, gradient, subspaces)
      !! `g_pq = 2 (F_qp - F_pq)`, zero on the redundant blocks
      !!
      !! The derivative of the energy with respect to `kappa_pq` under
      !! `C -> C exp(kappa)`. See the sign note in the module header.
      !!
      !! Antisymmetric by construction, which is what makes it a gradient with
      !! respect to an antisymmetric parametrisation: `kappa_pq` and
      !! `kappa_qp` are not independent, so the derivative cannot be either.
      type(mcscf_fock_t), intent(in) :: fock
      integer, intent(in) :: n_inactive, n_active
      real(dp), allocatable, intent(out) :: gradient(:, :)
      integer, intent(in), optional :: subspaces(:)
         !! Restricted-space partition; absent means every active-active
         !! rotation is redundant, which is the complete-space case

      integer :: n_mo, p, q

      n_mo = size(fock%general, 1)
      allocate (gradient(n_mo, n_mo))
      gradient = 0.0_dp
      do q = 1, n_mo
         do p = 1, n_mo
            if (is_redundant(p, q, n_inactive, n_active, subspaces)) cycle
            gradient(p, q) = 2.0_dp*(fock%general(q, p) - fock%general(p, q))
         end do
      end do
   end subroutine orbital_gradient

   subroutine rotation_parameters(n_mo, n_inactive, n_active, rows, cols, subspaces)
      !! The rotations that are real parameters, as a flat list
      !!
      !! Everything second order works in this list rather than in the `n_mo` by
      !! `n_mo` matrix the gradient arrives in: carry the redundant rotations
      !! along and the Hessian is singular by construction. Only `p > q`
      !! appears, since `kappa` is antisymmetric and the two triangles are the
      !! same variable seen twice.
      integer, intent(in) :: n_mo, n_inactive, n_active
      integer, allocatable, intent(out) :: rows(:), cols(:)
         !! `(rows(k), cols(k))` is the orbital pair parameter `k` rotates
      integer, intent(in), optional :: subspaces(:)   !! As `orbital_gradient`

      integer :: p, q, n_param

      n_param = 0
      do q = 1, n_mo
         do p = q + 1, n_mo
            if (.not. is_redundant(p, q, n_inactive, n_active, subspaces)) then
               n_param = n_param + 1
            end if
         end do
      end do

      allocate (rows(n_param), cols(n_param))
      n_param = 0
      do q = 1, n_mo
         do p = q + 1, n_mo
            if (is_redundant(p, q, n_inactive, n_active, subspaces)) cycle
            n_param = n_param + 1
            rows(n_param) = p
            cols(n_param) = q
         end do
      end do
   end subroutine rotation_parameters

   subroutine orbital_hessian_apply(this, x, hx, error)
      !! The orbital Hessian times one vector of rotation parameters
      class(orbital_hessian_operator_t), intent(inout) :: this
      real(dp), intent(in) :: x(:)      !! (n_param)
      real(dp), intent(out) :: hx(:)    !! (n_param)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: kappa_stack(:, :, :), fock_stack(:, :, :)
      integer :: n_mo, l

      if (error%has_error()) return
      n_mo = size(this%a_block, 1)
      allocate (kappa_stack(n_mo, n_mo, 1), fock_stack(n_mo, n_mo, 1))
      kappa_stack = 0.0_dp
      do l = 1, size(this%rows)
         kappa_stack(this%rows(l), this%cols(l), 1) = x(l)
         kappa_stack(this%cols(l), this%rows(l), 1) = -x(l)
      end do
      call one_index_fock_many(this%a_block, this%b_block, this%fock, this%dm1, this%dm2, &
                               this%n_inactive, this%n_active, kappa_stack, fock_stack)
      do l = 1, size(this%rows)
         hx(l) = 2.0_dp*(fock_stack(this%cols(l), this%rows(l), 1) &
                         - fock_stack(this%rows(l), this%cols(l), 1))
      end do
   end subroutine orbital_hessian_apply

   function approximate_hessian_diagonal(fock, rows, cols) result(diagonal)
      !! The one-electron approximation to the orbital Hessian's diagonal
      !!
      !!     H_pq,pq ~ 2 (n_q f_pp + n_p f_qq) - 2 (F_pp + F_qq)
      !!
      !! with `n` the occupations, `f = FI + FA` and `F` the generalised Fock.
      !! For an inactive-virtual pair it is `4 (f_aa - f_ii)`.
      type(mcscf_fock_t), intent(in) :: fock
      integer, intent(in) :: rows(:), cols(:)
      real(dp) :: diagonal(size(rows))

      integer :: l, p, q

      do l = 1, size(rows)
         p = rows(l)
         q = cols(l)
         diagonal(l) = 2.0_dp*(fock%occupation(q)*(fock%inactive(p, p) + fock%active(p, p)) &
                               + fock%occupation(p)*(fock%inactive(q, q) + fock%active(q, q))) &
                       - 2.0_dp*(fock%general(p, p) + fock%general(q, q))
      end do
   end function approximate_hessian_diagonal

   subroutine iterative_newton_step(operator, gradient, escape, kappa, lowest, predicted, &
                                    products, error)
      !! `newton_step` without building or diagonalising the Hessian
      !!
      !! `mqc_orbital_rotation`'s `subspace_newton_step` on the operator, with
      !! `approximate_hessian_diagonal` as its preconditioner, so the level
      !! shift and saddle escape are `level_shifted_step`'s on the projected
      !! Hessian. `lowest` is then an upper bound on the smallest curvature.
      type(orbital_hessian_operator_t), intent(inout) :: operator
      real(dp), intent(in) :: gradient(:, :)   !! (n_mo, n_mo), as `orbital_gradient`
      real(dp), intent(in) :: escape
         !! How far to displace a mode that has negative curvature and no
         !! gradient, in radians
      real(dp), allocatable, intent(out) :: kappa(:, :)
      real(dp), intent(out) :: lowest
      real(dp), intent(out) :: predicted
         !! What the quadratic model says the step is worth, as a positive
         !! energy decrease
      integer, intent(out) :: products
         !! Hessian-vector products spent
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat_gradient(:), step(:)
      integer :: n_mo, n_param, l

      lowest = 0.0_dp
      predicted = 0.0_dp
      products = 0
      if (error%has_error()) return
      n_mo = size(gradient, 1)
      n_param = size(operator%rows)
      allocate (kappa(n_mo, n_mo))
      kappa = 0.0_dp
      if (n_param == 0) return

      allocate (flat_gradient(n_param))
      do l = 1, n_param
         flat_gradient(l) = gradient(operator%rows(l), operator%cols(l))
      end do
      call subspace_newton_step(operator, &
                                approximate_hessian_diagonal(operator%fock, operator%rows, &
                                                             operator%cols), &
                                flat_gradient, escape, step, lowest, predicted, products, error, &
                                max_subspace=ITERATIVE_SUBSPACE, tolerance=ITERATIVE_TOLERANCE)
      if (error%has_error()) return

      do l = 1, n_param
         kappa(operator%rows(l), operator%cols(l)) = step(l)
         kappa(operator%cols(l), operator%rows(l)) = -step(l)
      end do
   end subroutine iterative_newton_step

   subroutine orbital_hessian(mol, orbitals, n_inactive, n_active, dm1, dm2, fock, &
                              rows, cols, hessian, error)
      !! The exact orbital Hessian at fixed CI, over the non-redundant rotations
      !!
      !! One column per parameter, each the differentiated Fock matrix of
      !! `one_index_fock` read through the gradient expression. Dense, because
      !! the point of building it is to *diagonalise* it: escaping a saddle
      !! needs the eigenvalues.
      !!
      !! The two integral blocks are `n_mo^2 n_occ^2` each. Everything an MCSCF
      !! Hessian needs has at most two general indices because the density
      !! matrices supply the rest, so no `n_mo^4` MO tensor is ever formed.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)      !! (n_ao, n_mo)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1(:, :), dm2(:, :, :, :)
      type(mcscf_fock_t), intent(in) :: fock
      integer, intent(in) :: rows(:), cols(:)     !! From `rotation_parameters`
      real(dp), allocatable, intent(out) :: hessian(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: a_mo(:, :, :, :), b_mo(:, :, :, :)

      if (error%has_error()) return
      if (size(rows) == 0) then
         allocate (hessian(0, 0))
         return
      end if
      call mo_integral_blocks(mol, orbitals, n_inactive + n_active, a_mo, b_mo)
      call orbital_hessian_from_blocks(a_mo, b_mo, n_inactive, n_active, dm1, dm2, fock, &
                                       rows, cols, hessian)
   end subroutine orbital_hessian

   subroutine orbital_hessian_from_blocks(a_block, b_block, n_inactive, n_active, dm1, dm2, &
                                          fock, rows, cols, hessian)
      !! `orbital_hessian`, from MO integral blocks the caller already holds
      real(dp), intent(in), contiguous :: a_block(:, :, :, :)
         !! (n_mo, n_occ, n_mo, n_occ): `(p q|r s)`, `q` and `s` occupied
      real(dp), intent(in), contiguous :: b_block(:, :, :, :)
         !! (n_mo, n_mo, n_occ, n_occ): `(p q|r s)`, `r` and `s` occupied
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1(:, :), dm2(:, :, :, :)
      type(mcscf_fock_t), intent(in) :: fock
      integer, intent(in) :: rows(:), cols(:)     !! From `rotation_parameters`
      real(dp), allocatable, intent(out) :: hessian(:, :)

      real(dp), allocatable :: kappa_block(:, :, :), transformed(:, :, :)
      integer :: n_mo, n_param, k, l, k0, k1, nk

      integer, parameter :: COLUMN_BLOCK = 48
         !! Columns per `one_index_fock_many` call: enough that its products
         !! are matrix-matrix, small enough that one thread's stack of
         !! `(n_mo, n_mo)` matrices stays a few megabytes.

      n_mo = size(a_block, 1)
      n_param = size(rows)
      allocate (hessian(n_param, n_param))
      if (n_param == 0) return

      !$omp parallel do schedule(dynamic) default(shared) &
      !$omp    private(k0, k1, nk, k, l, kappa_block, transformed)
      do k0 = 1, n_param, COLUMN_BLOCK
         k1 = min(k0 + COLUMN_BLOCK - 1, n_param)
         nk = k1 - k0 + 1
         allocate (kappa_block(n_mo, n_mo, nk), transformed(n_mo, n_mo, nk))
         kappa_block = 0.0_dp
         do k = k0, k1
            kappa_block(rows(k), cols(k), k - k0 + 1) = 1.0_dp
            kappa_block(cols(k), rows(k), k - k0 + 1) = -1.0_dp
         end do
         call one_index_fock_many(a_block, b_block, fock, dm1, dm2, n_inactive, n_active, &
                                  kappa_block, transformed)
         do k = k0, k1
            do l = 1, n_param
               hessian(l, k) = 2.0_dp*(transformed(cols(l), rows(l), k - k0 + 1) &
                                       - transformed(rows(l), cols(l), k - k0 + 1))
            end do
         end do
         deallocate (kappa_block, transformed)
      end do
      !$omp end parallel do

      ! Symmetric in exact arithmetic. Away from a stationary point the two
      ! triangles differ by rounding and by the antisymmetric term that
      ! distinguishes differentiating the gradient from differentiating the
      ! energy twice, which this average removes.
      hessian = 0.5_dp*(hessian + transpose(hessian))
   end subroutine orbital_hessian_from_blocks

   subroutine mo_integral_blocks(mol, orbitals, n_occ, a_block, b_block)
      !! The two MO integral blocks the orbital Hessian and the generalised
      !! Fock matrix are built from, from one pass over the AO integrals
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)      !! (n_ao, n_mo)
      integer, intent(in) :: n_occ
      real(dp), allocatable, intent(out) :: a_block(:, :, :, :)
         !! (n_mo, n_occ, n_mo, n_occ): `(p q|r s)`, `q` and `s` occupied
      real(dp), allocatable, intent(out) :: b_block(:, :, :, :)
         !! (n_mo, n_mo, n_occ, n_occ): `(p q|r s)`, `r` and `s` occupied

      real(dp), allocatable :: eri_packed(:, :)

      call mol%eris_packed(eri_packed)
      call transform_block(eri_packed, orbitals, orbitals(:, 1:n_occ), orbitals, &
                           orbitals(:, 1:n_occ), a_block)
      call transform_block(eri_packed, orbitals, orbitals, orbitals(:, 1:n_occ), &
                           orbitals(:, 1:n_occ), b_block)
   end subroutine mo_integral_blocks

   subroutine fock_from_blocks(mol, orbitals, n_inactive, n_active, dm1, dm2, a_block, &
                               b_block, fock)
      !! `generalized_fock`, from the MO integral blocks instead of direct AO
      !! Fock builds
      !!
      !! `FI = h + sum_i [2 (pq|ii) - (pi|qi)]` and
      !! `FA = sum_tu D_tu [(pq|tu) - (pt|qu)/2]` over every orbital pair, and
      !! `(n u|v w)` is a slice of `a_block`.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)      !! (n_ao, n_mo)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1(:, :)           !! Active one-particle density
      real(dp), intent(in) :: dm2(:, :, :, :)     !! Active two-particle density
      real(dp), intent(in) :: a_block(:, :, :, :)
         !! (n_mo, n_occ, n_mo, n_occ): `(p q|r s)`, `q` and `s` occupied
      real(dp), intent(in) :: b_block(:, :, :, :)
         !! (n_mo, n_mo, n_occ, n_occ): `(p q|r s)`, `r` and `s` occupied
      type(mcscf_fock_t), intent(out) :: fock

      real(dp), allocatable :: h_ao(:, :), work(:, :)
      real(dp) :: accumulated
      integer :: n_ao, n_mo, n_occ, i, t, u, v, w, n, ta, ua

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      n_occ = n_inactive + n_active

      call mol%core_hamiltonian(h_ao)
      allocate (fock%inactive(n_mo, n_mo), fock%active(n_mo, n_mo), work(n_ao, n_mo))
      call pic_gemm(h_ao, orbitals, work)
      call pic_gemm(orbitals, work, fock%inactive, transa="T")
      do i = 1, n_inactive
         fock%inactive = fock%inactive + 2.0_dp*b_block(:, :, i, i) - a_block(:, i, :, i)
      end do

      fock%active = 0.0_dp
      do u = 1, n_active
         ua = n_inactive + u
         do t = 1, n_active
            ta = n_inactive + t
            fock%active = fock%active + dm1(t, u)*(b_block(:, :, ta, ua) &
                                                   - 0.5_dp*a_block(:, ta, :, ua))
         end do
      end do

      allocate (fock%general(n_mo, n_mo), fock%occupation(n_mo))
      fock%general = 0.0_dp
      fock%occupation = 0.0_dp

      do i = 1, n_inactive
         fock%occupation(i) = 2.0_dp
         do n = 1, n_mo
            fock%general(i, n) = 2.0_dp*(fock%inactive(n, i) + fock%active(n, i))
         end do
      end do

      do t = 1, n_active
         fock%occupation(n_inactive + t) = dm1(t, t)
         do n = 1, n_mo
            accumulated = 0.0_dp
            do u = 1, n_active
               accumulated = accumulated + dm1(t, u)*fock%inactive(n, n_inactive + u)
            end do
            do w = 1, n_active
               do v = 1, n_active
                  do u = 1, n_active
                     accumulated = accumulated + dm2(t, u, v, w)* &
                                   a_block(n, n_inactive + u, n_inactive + v, n_inactive + w)
                  end do
               end do
            end do
            fock%general(n_inactive + t, n) = accumulated
         end do
      end do
   end subroutine fock_from_blocks

   subroutine one_index_fock(a_block, b_block, fock, dm1, dm2, n_inactive, n_active, &
                             kappa, transformed)
      !! The generalised Fock matrix differentiated along one orbital rotation
      !!
      !!     (H kappa)_pq = 2 (Ft_qp - Ft_pq)
      !!
      !! with `Ft` built exactly as `generalized_fock` builds `F`, but with each
      !! orbital in turn replaced by `C kappa`. `one_index_fock_many` for one
      !! rotation.
      !!
      !! The density matrices are the ones the CI produced and are held fixed:
      !! this is the orbital-orbital block at fixed CI, not the full second
      !! derivative of the two-step energy. The difference is a Schur complement
      !! that only makes eigenvalues smaller, so a negative direction found here
      !! is a genuine one.
      real(dp), intent(in), contiguous :: a_block(:, :, :, :)
         !! `(p q|r s)`, `q` and `s` occupied
      real(dp), intent(in), contiguous :: b_block(:, :, :, :)
         !! `(p q|r s)`, `r` and `s` occupied
      type(mcscf_fock_t), intent(in) :: fock
      real(dp), intent(in) :: dm1(:, :), dm2(:, :, :, :)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: kappa(:, :)           !! (n_mo, n_mo), antisymmetric
      real(dp), intent(out) :: transformed(:, :)    !! (n_mo, n_mo)

      real(dp), allocatable :: kappa_stack(:, :, :), fock_stack(:, :, :)
      integer :: n_mo

      n_mo = size(kappa, 1)
      allocate (kappa_stack(n_mo, n_mo, 1), fock_stack(n_mo, n_mo, 1))
      kappa_stack(:, :, 1) = kappa
      call one_index_fock_many(a_block, b_block, fock, dm1, dm2, n_inactive, n_active, &
                               kappa_stack, fock_stack)
      transformed = fock_stack(:, :, 1)
   end subroutine one_index_fock

   subroutine one_index_fock_many(a_block, b_block, fock, dm1, dm2, n_inactive, n_active, &
                                  kappas, transformed)
      !! `one_index_fock` for a stack of rotations, with the integral
      !! contractions done as matrix products over the whole stack
      !!
      !! Uses that each `kappa` is antisymmetric and that the Fock and
      !! occupation matrices are symmetric, so `kappa M = -(M kappa)^T`.
      real(dp), intent(in), contiguous :: a_block(:, :, :, :)
         !! `(p q|r s)`, `q` and `s` occupied
      real(dp), intent(in), contiguous :: b_block(:, :, :, :)
         !! `(p q|r s)`, `r` and `s` occupied
      type(mcscf_fock_t), intent(in) :: fock
      real(dp), intent(in) :: dm1(:, :), dm2(:, :, :, :)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in), contiguous :: kappas(:, :, :)
         !! (n_mo, n_mo, n_set), each antisymmetric
      real(dp), intent(out), contiguous :: transformed(:, :, :)    !! (n_mo, n_mo, n_set)

      call fock_kernel(size(kappas, 1), n_inactive, n_active, size(kappas, 3), a_block, &
                       b_block, fock, dm1, dm2, kappas, transformed)
   end subroutine one_index_fock_many

   subroutine fock_kernel(n_mo, n_inactive, n_active, n_set, a_block, b_block, fock, dm1, &
                          dm2, kappas, transformed)
      !! `one_index_fock_many`, with explicit shapes so that slices of the
      !! integral blocks go to BLAS as an element and a leading dimension
      integer, intent(in) :: n_mo, n_inactive, n_active, n_set
      real(dp), intent(in) :: a_block(n_mo, n_inactive + n_active, n_mo, n_inactive + n_active)
      real(dp), intent(in) :: b_block(n_mo, n_mo, n_inactive + n_active, n_inactive + n_active)
      type(mcscf_fock_t), intent(in) :: fock
      real(dp), intent(in) :: dm1(n_active, n_active)
      real(dp), intent(in) :: dm2(n_active, n_active, n_active, n_active)
      real(dp), intent(in) :: kappas(n_mo, n_mo, n_set)
      real(dp), intent(out) :: transformed(n_mo, n_mo, n_set)

      real(dp), allocatable :: occupations(:, :), product(:, :)
      real(dp), allocatable :: rotated_densities(:, :, :), rotated_potentials(:, :, :)
      real(dp), allocatable :: fi_kappa(:, :), fa_kappa(:, :)
      real(dp), allocatable :: f_inactive(:, :), f_active(:, :)
      real(dp), allocatable :: a_gaaa(:, :), eri_gaaa(:, :, :, :, :)
      real(dp), allocatable :: active_rows(:, :), dm2_rows(:, :)
      integer :: n_occ, n_act3, k, i, u, v, w, lo, hi, ia

      n_occ = n_inactive + n_active
      n_act3 = n_active**3
      ia = n_inactive + 1

      ! The inactive density is C_i C_i^T, so its derivative replaces one factor
      ! at a time and is non-zero only where exactly one index is inactive --
      ! which is also why a rotation between two inactive orbitals moves
      ! nothing. The active one is `kappa D - D kappa`.
      allocate (occupations(n_mo, n_mo), product(n_mo, n_mo*n_set))
      occupations = 0.0_dp
      occupations(n_inactive + 1:n_occ, n_inactive + 1:n_occ) = dm1
      call pic_dgemm_x("N", "N", n_mo, n_mo*n_set, n_mo, 1.0_dp, occupations, n_mo, &
                       kappas, n_mo, 0.0_dp, product, n_mo)
      allocate (rotated_densities(n_mo, n_mo, 2*n_set))
      do k = 1, n_set
         lo = (k - 1)*n_mo + 1
         hi = k*n_mo
         rotated_densities(:, :, 2*k - 1) = 0.0_dp
         rotated_densities(:, 1:n_inactive, 2*k - 1) = 2.0_dp*kappas(:, 1:n_inactive, k)
         rotated_densities(1:n_inactive, :, 2*k - 1) = rotated_densities(1:n_inactive, :, 2*k - 1) &
                                                       - 2.0_dp*kappas(1:n_inactive, :, k)
         rotated_densities(:, :, 2*k) = -transpose(product(:, lo:hi)) - product(:, lo:hi)
      end do
      deallocate (product)

      allocate (rotated_potentials(n_mo, n_occ, 2*n_set))
      call potential_kernel(n_mo, n_occ, 2*n_set, a_block, b_block, rotated_densities, rotated_potentials)
      deallocate (rotated_densities)

      ! Two contributions to each Fock matrix: the orbitals it is expressed in,
      ! which give a commutator, and the orbitals its density was built from,
      ! which give the potential above.
      allocate (fi_kappa(n_mo, n_mo*n_set), fa_kappa(n_mo, n_mo*n_set))
      call pic_dgemm_x("N", "N", n_mo, n_mo*n_set, n_mo, 1.0_dp, fock%inactive, n_mo, &
                       kappas, n_mo, 0.0_dp, fi_kappa, n_mo)
      call pic_dgemm_x("N", "N", n_mo, n_mo*n_set, n_mo, 1.0_dp, fock%active, n_mo, &
                       kappas, n_mo, 0.0_dp, fa_kappa, n_mo)

      ! (nu|vw) has four orbitals in it and so four terms, one per index, each
      ! a product of `kappa` with a slice of an integral block written straight
      ! into a strided view of `eri_gaaa(n, u, v, w, k)`.
      allocate (a_gaaa(n_mo, n_act3), eri_gaaa(n_mo, n_active, n_active, n_active, n_set))
      a_gaaa = reshape(a_block(:, ia:n_occ, ia:n_occ, ia:n_occ), [n_mo, n_act3])
      do k = 1, n_set
         ! sum_m kappa_mn (m u|v w): the first index transformed
         call pic_dgemm_x("T", "N", n_mo, n_act3, n_mo, 1.0_dp, kappas(1, 1, k), n_mo, &
                          a_gaaa, n_mo, 0.0_dp, eri_gaaa(1, 1, 1, 1, k), n_mo)
         do w = 1, n_active
            do v = 1, n_active
               ! sum_m (n m|v w) kappa_mu: the second index
               call pic_dgemm_x("N", "N", n_mo, n_active, n_mo, 1.0_dp, &
                                b_block(1, 1, n_inactive + v, n_inactive + w), n_mo, &
                                kappas(1, ia, k), n_mo, 1.0_dp, eri_gaaa(1, 1, v, w, k), n_mo)
            end do
         end do
         do w = 1, n_active
            do u = 1, n_active
               ! sum_m (n u|m w) kappa_mv: the third index
               call pic_dgemm_x("N", "N", n_mo, n_active, n_mo, 1.0_dp, &
                                a_block(1, n_inactive + u, 1, n_inactive + w), n_mo*n_occ, &
                                kappas(1, ia, k), n_mo, 1.0_dp, eri_gaaa(1, u, 1, w, k), &
                                n_mo*n_active)
            end do
         end do
         do v = 1, n_active
            do u = 1, n_active
               ! sum_m (n u|v m) kappa_mw = sum_m (n u|m v) kappa_mw: the fourth
               call pic_dgemm_x("N", "N", n_mo, n_active, n_mo, 1.0_dp, &
                                a_block(1, n_inactive + u, 1, n_inactive + v), n_mo*n_occ, &
                                kappas(1, ia, k), n_mo, 1.0_dp, eri_gaaa(1, u, v, 1, k), &
                                n_mo*n_active*n_active)
            end do
         end do
      end do
      deallocate (a_gaaa)

      ! The rows, assembled exactly as `generalized_fock` assembles them. The
      ! virtual rows stay zero because the density is, and the density is not
      ! what is being differentiated.
      allocate (f_inactive(n_mo, n_occ), f_active(n_mo, n_occ))
      allocate (active_rows(n_mo, n_active), dm2_rows(n_active, n_act3))
      dm2_rows = reshape(dm2, [n_active, n_act3])
      do k = 1, n_set
         lo = (k - 1)*n_mo + 1
         hi = k*n_mo
         f_inactive = fi_kappa(:, lo:lo + n_occ - 1) &
                      + transpose(fi_kappa(1:n_occ, lo:hi)) + rotated_potentials(:, :, 2*k - 1)
         f_active = fa_kappa(:, lo:lo + n_occ - 1) &
                    + transpose(fa_kappa(1:n_occ, lo:hi)) + rotated_potentials(:, :, 2*k)

         transformed(:, :, k) = 0.0_dp
         do i = 1, n_inactive
            transformed(i, :, k) = 2.0_dp*(f_inactive(:, i) + f_active(:, i))
         end do
         if (n_active > 0) then
            call pic_gemm(f_inactive(:, n_inactive + 1:n_occ), dm1, active_rows, transb="T")
            call pic_dgemm_x("N", "T", n_mo, n_active, n_act3, 1.0_dp, eri_gaaa(1, 1, 1, 1, k), &
                             n_mo, dm2_rows, n_active, 1.0_dp, active_rows, n_mo)
            transformed(n_inactive + 1:n_occ, :, k) = transpose(active_rows)
         end if
      end do

      deallocate (occupations, rotated_potentials, fi_kappa, fa_kappa, f_inactive, f_active)
      deallocate (eri_gaaa, active_rows, dm2_rows)
   end subroutine fock_kernel

   subroutine transformed_potential(a_block, b_block, n_occ, density, potential)
      !! `J - K/2` of a differentiated density, in the MO basis
      !!
      !! Only the occupied columns are produced, because those are the only ones
      !! the generalised Fock reads: its virtual rows are zero, so nothing ever
      !! asks what the potential does between two empty orbitals.
      !!
      !! **The density must vanish on the virtual-virtual block**, which is what
      !! makes the two integral blocks sufficient. Both densities this is used
      !! for are one-index transforms of a density carried by occupied orbitals,
      !! so one index is always occupied. A density without that structure would
      !! need `(virtual virtual|virtual virtual)` integrals, which are not here.
      real(dp), intent(in), contiguous :: a_block(:, :, :, :)
         !! `(p q|r s)`, `q` and `s` occupied
      real(dp), intent(in), contiguous :: b_block(:, :, :, :)
         !! `(p q|r s)`, `r` and `s` occupied
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: density(:, :)         !! (n_mo, n_mo), symmetric
      real(dp), intent(out) :: potential(:, :)      !! (n_mo, n_occ)

      real(dp), allocatable :: density_stack(:, :, :), potential_stack(:, :, :)
      integer :: n_mo

      n_mo = size(density, 1)
      allocate (density_stack(n_mo, n_mo, 1), potential_stack(n_mo, n_occ, 1))
      density_stack(:, :, 1) = density
      call transformed_potential_many(a_block, b_block, n_occ, density_stack, potential_stack)
      potential = potential_stack(:, :, 1)
   end subroutine transformed_potential

   subroutine transformed_potential_many(a_block, b_block, n_occ, densities, potentials)
      !! `transformed_potential` for a stack of densities, as matrix products
      !! over the whole stack
      real(dp), intent(in), contiguous :: a_block(:, :, :, :)
         !! (n_mo, n_occ, n_mo, n_occ): `(p q|r s)`, `q` and `s` occupied
      real(dp), intent(in), contiguous :: b_block(:, :, :, :)
         !! (n_mo, n_mo, n_occ, n_occ): `(p q|r s)`, `r` and `s` occupied
      integer, intent(in) :: n_occ
      real(dp), intent(in), contiguous :: densities(:, :, :)
         !! (n_mo, n_mo, n_set), each symmetric
      real(dp), intent(out), contiguous :: potentials(:, :, :)   !! (n_mo, n_occ, n_set)

      call potential_kernel(size(densities, 1), n_occ, size(densities, 3), a_block, b_block, &
                            densities, potentials)
   end subroutine transformed_potential_many

   subroutine potential_kernel(n_mo, n_occ, n_set, a_block, b_block, densities, potentials)
      !! `transformed_potential_many`, with explicit shapes so that slices of
      !! the integral blocks go to BLAS as an element and a leading dimension
      integer, intent(in) :: n_mo, n_occ, n_set
      real(dp), intent(in) :: a_block(n_mo, n_occ, n_mo, n_occ)
      real(dp), intent(in) :: b_block(n_mo, n_mo, n_occ, n_occ)
      real(dp), intent(in) :: densities(n_mo, n_mo, n_set)
      real(dp), intent(out) :: potentials(n_mo, n_occ, n_set)

      real(dp), allocatable :: occupied(:, :), weight(:, :), coulomb(:, :)
      real(dp), allocatable :: exchange(:, :, :), virtual(:, :, :), product(:, :)
      real(dp), allocatable :: mine(:, :, :)
      integer :: m, n_virt, k, q, r, s, p, row0, rows

      integer, parameter :: ROW_BLOCK = 512
         !! Rows of the Coulomb product per thread chunk

      m = n_mo*n_occ
      n_virt = n_mo - n_occ

      ! Coulomb. As a matrix, `a_block` is `(p q|r s)` with rows `(p q)` and
      ! columns `(r s)`. The block has only the second density index occupied;
      ! a pair with the first index occupied is the same integral read the
      ! other way round, so doubling the rows that are not occupied covers both.
      allocate (occupied(m, n_set), weight(m, n_set), coulomb(m, n_set))
      do k = 1, n_set
         occupied(:, k) = reshape(densities(:, 1:n_occ, k), [m])
         weight(:, k) = occupied(:, k)
         do s = 1, n_occ
            weight((s - 1)*n_mo + n_occ + 1:s*n_mo, k) = &
               2.0_dp*weight((s - 1)*n_mo + n_occ + 1:s*n_mo, k)
         end do
      end do
      ! Threaded over row blocks with sequential BLAS. Called from inside the
      ! orbital Hessian's own parallel loop, these regions run on one thread.
      !$omp parallel do schedule(static) default(shared) private(row0, rows)
      do row0 = 1, m, ROW_BLOCK
         rows = min(ROW_BLOCK, m - row0 + 1)
         call pic_dgemm_x("N", "N", rows, n_set, m, 1.0_dp, a_block(row0, 1, 1, 1), m, &
                          weight, m, 0.0_dp, coulomb(row0, 1), m)
      end do
      !$omp end parallel do

      ! Exchange, the half with both density indices in the occupied columns:
      ! `b_block(:, :, :, q)` is `(p r|s q)` with columns `(r s)`.
      allocate (exchange(n_mo, n_occ, n_set))
      !$omp parallel do schedule(static) default(shared)
      do q = 1, n_occ
         call pic_dgemm_x("N", "N", n_mo, n_set, m, 1.0_dp, b_block(1, 1, 1, q), n_mo, &
                          occupied, m, 0.0_dp, exchange(1, q, 1), n_mo*n_occ)
      end do
      !$omp end parallel do

      ! Exchange, the half with a virtual second density index. Swapping the
      ! density indices moves them to different slots of the integral, so it
      ! comes from the other block: sum over occupied r and virtual s of
      ! D_rs (s q|p r). The virtual rows of `a_block(:, :, :, r)` are a matrix
      ! over `(q p)` with leading dimension `n_mo`.
      if (n_virt > 0) then
         allocate (virtual(n_virt, n_set, n_occ))
         do r = 1, n_occ
            do k = 1, n_set
               virtual(:, k, r) = densities(r, n_occ + 1:n_mo, k)
            end do
         end do
         !$omp parallel default(shared) private(r, k, p, product, mine)
         allocate (product(m, n_set), mine(n_mo, n_occ, n_set))
         mine = 0.0_dp
         !$omp do schedule(static)
         do r = 1, n_occ
            call pic_dgemm_x("T", "N", m, n_set, n_virt, 1.0_dp, a_block(n_occ + 1, 1, 1, r), &
                             n_mo, virtual(1, 1, r), n_virt, 0.0_dp, product, m)
            do k = 1, n_set
               do p = 1, n_mo
                  mine(p, :, k) = mine(p, :, k) + product((p - 1)*n_occ + 1:p*n_occ, k)
               end do
            end do
         end do
         !$omp end do
         !$omp critical
         exchange = exchange + mine
         !$omp end critical
         deallocate (product, mine)
         !$omp end parallel
         deallocate (virtual)
      end if

      do k = 1, n_set
         potentials(:, :, k) = reshape(coulomb(:, k), [n_mo, n_occ]) - 0.5_dp*exchange(:, :, k)
      end do
      deallocate (occupied, weight, coulomb, exchange)
   end subroutine potential_kernel

   subroutine newton_step(hessian, gradient, rows, cols, escape, kappa, lowest, &
                          predicted, error)
      !! The Newton step for an MCSCF, as an antisymmetric `kappa`
      !!
      !! The step itself is `mqc_orbital_rotation`'s `level_shifted_step` --
      !! level shifting, the saddle escape and the predicted gain all live
      !! there, shared with the second-order SCF. What is left here is the only
      !! part that is MCSCF's: `(rows, cols)` says which orbital pair each
      !! parameter rotates, so the gradient is gathered out of the `n_mo` by
      !! `n_mo` matrix `orbital_gradient` produces and the step scattered back
      !! into an antisymmetric matrix `rotation_matrix` can exponentiate.
      real(dp), intent(in) :: hessian(:, :)
      real(dp), intent(in) :: gradient(:, :)   !! (n_mo, n_mo), as `orbital_gradient`
      integer, intent(in) :: rows(:), cols(:)  !! From `rotation_parameters`
      real(dp), intent(in) :: escape
         !! How far to displace a mode that has negative curvature and no
         !! gradient, in radians
      real(dp), allocatable, intent(out) :: kappa(:, :)
      real(dp), intent(out) :: lowest
         !! Smallest Hessian eigenvalue. Negative means the point the step was
         !! taken from is not a minimum, whatever the gradient says.
      real(dp), intent(out) :: predicted
         !! What the quadratic model says the step is worth, as a positive
         !! energy decrease.
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: flat_gradient(:), step(:)
      integer :: n_mo, n_param, l

      lowest = 0.0_dp
      predicted = 0.0_dp
      if (error%has_error()) return
      n_mo = size(gradient, 1)
      n_param = size(rows)
      allocate (kappa(n_mo, n_mo))
      kappa = 0.0_dp
      if (n_param == 0) return

      allocate (flat_gradient(n_param))
      do l = 1, n_param
         flat_gradient(l) = gradient(rows(l), cols(l))
      end do

      call level_shifted_step(hessian, flat_gradient, escape, step, lowest, &
                              predicted, error)
      if (error%has_error()) return

      do l = 1, n_param
         kappa(rows(l), cols(l)) = step(l)
         kappa(cols(l), rows(l)) = -step(l)
      end do

      deallocate (flat_gradient, step)
   end subroutine newton_step

   subroutine run_czt_casscf(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                             result, error, max_iterations, gradient_tol, verbose, &
                             subspaces, min_electrons, max_electrons, n_states, weights, &
                             iterative_hessian)
      !! Two-step CASSCF: solve the CI, move the orbitals, repeat
      !!
      !! Each macro-iteration solves the CI problem exactly at the current
      !! orbitals, builds the density matrices from it, and takes one Newton
      !! step on the orbitals. The CI is re-solved from the previous vector,
      !! which after the first few iterations is nearly the answer already.
      !!
      !! Two-step, so the orbital step ignores the coupling between orbital
      !! rotations and the CI coefficients. That costs iterations near the
      !! solution and nothing else, because the neglected block only makes the
      !! true curvature smaller than the one used here.
      !!
      !! **Convergence means a small gradient *and* no negative curvature.** An
      !! optimiser built from the gradient alone stops wherever it vanishes,
      !! which on a symmetric starting guess can be a saddle: the rotations that
      !! would break the symmetry have exactly zero gradient and keep it. So the
      !! Hessian is built and diagonalised even on the iteration that looks
      !! converged.
      !!
      !! There is no extrapolation. Against a Newton step DIIS was measured
      !! worth one iteration either way, for an extra CI solve per iteration.
      ! TODO(mqc): `max_iterations = 0` skips the loop, leaving `largest`
      ! undefined and `dm1`/`dm2` unallocated where the assignments below read
      ! all three.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      type(casscf_result_t), intent(out) :: result
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: gradient_tol
      logical, intent(in), optional :: verbose
      integer, intent(in), optional :: subspaces(:), min_electrons(:), max_electrons(:)
         !! An occupation-restricted space, as `keywords.mcscf.ormas` gives it.
         !! Absent is a complete active space. All three or none.
      integer, intent(in), optional :: n_states
         !! State-averaged CASSCF over this many singlet roots when greater
         !! than 1: the orbitals minimise `E_SA = sum_J w_J E_J`, and the
         !! orbital gradient and Hessian are built from the weighted densities.
         !! Needs `weights`, `n_alpha == n_beta`, and a complete active space.
      real(dp), intent(in), optional :: weights(:)
         !! One weight per state, summing to one
      logical, intent(in), optional :: iterative_hessian
         !! Take each Newton step iteratively (`iterative_newton_step`) instead
         !! of building and diagonalising the orbital Hessian. Absent, it is
         !! iterative above `ITERATIVE_HESSIAN_ABOVE` rotations.

      type(ormas_space_t) :: space
      logical :: restricted
      type(casci_result_t) :: ci, trial_ci
      type(mcscf_fock_t) :: fock
      type(link_table_t) :: alpha, beta
      real(dp), allocatable, target :: a_mo(:, :, :, :), b_mo(:, :, :, :)
      integer :: hessian_products
      type(orbital_hessian_operator_t) :: hessian_operator
      logical :: iterative
      real(dp), allocatable :: current(:, :), updated(:, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp), allocatable :: gradient(:, :), hessian(:, :), kappa(:, :), rotation(:, :)
      real(dp), allocatable :: step_matrix(:, :)
      real(dp), allocatable :: guess(:, :, :)
      real(dp), allocatable :: flat_guess(:, :)
      integer, allocatable :: rows(:), cols(:)
      character(len=160) :: line
      real(dp) :: tol, largest, previous, trust, scaling, lowest, predicted
      real(dp) :: energy, trial_energy
      integer :: cycles, iteration, n_ao, n_mo, trial, n_roots, j
      logical :: loud, have_guess, accepted, saddle, averaged

      integer, parameter :: MAX_BACKTRACKS = 12

      if (error%has_error()) return

      n_roots = 1
      if (present(n_states)) n_roots = n_states
      averaged = n_roots > 1
      if (averaged) then
         if (.not. present(weights)) then
            call error%set(ERROR_VALIDATION, "mcscf: n_states > 1 needs weights")
            return
         end if
         if (size(weights) /= n_roots) then
            call error%set(ERROR_VALIDATION, "mcscf: "//to_char(size(weights))// &
                           " weights for "//to_char(n_roots)//" states")
            return
         end if
         if (n_alpha /= n_beta) then
            call error%set(ERROR_VALIDATION, "mcscf: state averaging needs a singlet "// &
                           "(equal active alpha and beta electrons); this active "// &
                           "space has "//to_char(n_alpha)//" alpha and "// &
                           to_char(n_beta)//" beta.")
            return
         end if
         if (present(subspaces)) then
            call error%set(ERROR_VALIDATION, "mcscf: state averaging over an "// &
                           "occupation-restricted (ORMAS) space is not implemented")
            return
         end if
      end if

      cycles = 50
      if (present(max_iterations)) cycles = max_iterations
      ! `keywords.mcscf.max_macro_iter` reaches here unchecked -- the schema
      ! allow-lists the key without a range -- and a zero skips the macro loop
      ! entirely, leaving `largest` undefined and `dm1`/`dm2` unallocated for
      ! the assembly below to read. Refusing beats reporting a gradient norm
      ! that was never computed.
      if (cycles < 1) then
         call error%set(ERROR_VALIDATION, "mcscf: max_macro_iter must be at least 1; "// &
                        "a run with no macro-iterations has no orbitals to report")
         return
      end if
      tol = 1.0e-6_dp
      if (present(gradient_tol)) tol = gradient_tol
      loud = .false.
      if (present(verbose)) loud = verbose

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      allocate (current(n_ao, n_mo), updated(n_ao, n_mo), step_matrix(n_mo, n_mo))
      current = orbitals

      restricted = present(subspaces)
      if (restricted) then
         call build_ormas_space(subspaces, n_active, n_alpha, n_beta, min_electrons, &
                                max_electrons, space, error)
         if (error%has_error()) return
      end if

      call build_link_table(n_active, n_alpha, alpha, error)
      call build_link_table(n_active, n_beta, beta, error)
      if (error%has_error()) return

      ! Which rotations are real parameters does not change as the orbitals
      ! move, so the list is built once.
      if (restricted) then
         call rotation_parameters(n_mo, n_inactive, n_active, rows, cols, subspaces)
      else
         call rotation_parameters(n_mo, n_inactive, n_active, rows, cols)
      end if

      iterative = size(rows) > ITERATIVE_HESSIAN_ABOVE
      if (present(iterative_hessian)) iterative = iterative_hessian
      if (iterative) then
         hessian_operator%n_inactive = n_inactive
         hessian_operator%n_active = n_active
         hessian_operator%rows = rows
         hessian_operator%cols = cols
      end if

      if (loud) then
         call logger%info("")
         if (restricted) then
            call logger%info("  occupation-restricted active space SCF")
         else if (averaged) then
            call logger%info("  state-averaged complete active space SCF, "// &
                             to_char(n_roots)//" states")
         else
            call logger%info("  complete active space SCF")
         end if
         call logger%info("    iter            energy          change      max gradient"// &
                          "       trust    curvature")
      end if

      have_guess = .false.
      trust = MAX_ROTATION

      ! The starting point, so the first step has something to improve on.
      if (restricted) then
         call solve_ci(mol, current, n_inactive, n_active, n_alpha, n_beta, alpha, beta, &
                       guess, have_guess, ci, error, subspaces, min_electrons, &
                       max_electrons, flat_guess)
      else
         call solve_ci(mol, current, n_inactive, n_active, n_alpha, n_beta, alpha, beta, &
                       guess, have_guess, ci, error, n_roots=n_roots)
      end if
      if (error%has_error()) return
      energy = ci%energy
      if (averaged) energy = weighted_sum(ci%energies, weights)
      previous = energy

      do iteration = 1, cycles
         result%iterations = iteration

         if (allocated(dm1)) deallocate (dm1)
         if (allocated(dm2)) deallocate (dm2)
         if (restricted) then
            call ormas_density_matrices(space, ci%ci_flat, dm1, dm2, error)
         else if (averaged) then
            call sa_density_matrices(ci%vectors, weights, alpha, beta, dm1, dm2, error)
         else
            call active_space_rdms(ci%ci_vector, alpha, beta, dm1, dm2, error)
         end if
         if (error%has_error()) return

         ! One AO integral pass per macro-iteration: the generalised Fock
         ! matrix and the Hessian are both read out of these blocks.
         call mo_integral_blocks(mol, current, n_inactive + n_active, a_mo, b_mo)
         call fock_from_blocks(mol, current, n_inactive, n_active, dm1, dm2, a_mo, b_mo, fock)
         if (allocated(gradient)) deallocate (gradient)
         if (restricted) then
            call orbital_gradient(fock, n_inactive, n_active, gradient, subspaces)
         else
            call orbital_gradient(fock, n_inactive, n_active, gradient)
         end if
         largest = maxval(abs(gradient))

         ! Built before the convergence test rather than after it, because the
         ! test needs the smallest eigenvalue: a zero gradient at a saddle is
         ! still a zero gradient.
         if (allocated(hessian)) deallocate (hessian)
         if (allocated(kappa)) deallocate (kappa)
         if (iterative) then
            hessian_operator%a_block => a_mo
            hessian_operator%b_block => b_mo
            hessian_operator%fock = fock
            hessian_operator%dm1 = dm1
            hessian_operator%dm2 = dm2
            call iterative_newton_step(hessian_operator, gradient, MAX_ROTATION, kappa, lowest, &
                                       predicted, hessian_products, error)
            nullify (hessian_operator%a_block, hessian_operator%b_block)
            deallocate (a_mo, b_mo)
         else
            call orbital_hessian_from_blocks(a_mo, b_mo, n_inactive, n_active, dm1, dm2, &
                                             fock, rows, cols, hessian)
            deallocate (a_mo, b_mo)
            call newton_step(hessian, gradient, rows, cols, MAX_ROTATION, kappa, &
                             lowest, predicted, error)
         end if
         if (error%has_error()) return

         if (loud) then
            write (line, "(a,i4,f20.12,2es16.4,2es12.2)") "    ", iteration, &
               energy, energy - previous, largest, trust, lowest
            call logger%info(trim(line))
         end if
         previous = energy

         saddle = lowest < SADDLE_CURVATURE
         if (largest < tol) then
            if (.not. saddle) then
               result%converged = .true.
               exit
            end if
            ! Leaving a saddle is a fresh direction, so it gets a fresh trust
            ! radius: the old one records how the previous direction behaved and
            ! is usually small by the time a run has settled.
            trust = MAX_ROTATION
            if (loud) call logger%info("    the gradient has vanished at a saddle "// &
                                       "point; following the negative curvature out")
         end if

         ! The trust region is what makes this converge rather than oscillate:
         ! an exact Hessian still describes the surface only near the point it
         ! was built at, and early steps are not near it. It is also what picks
         ! the sign of a saddle escape, where the second-order model is
         ! indifferent between the two directions.
         accepted = .false.

         ! ---- take it, backtracking until it descends ----------------------
         do trial = 1, MAX_BACKTRACKS
            if (accepted) exit
            scaling = maxval(abs(kappa))
            if (scaling > trust) then
               step_matrix = kappa*(trust/scaling)
            else
               step_matrix = kappa
            end if

            if (allocated(rotation)) deallocate (rotation)
            call rotation_matrix(step_matrix, rotation)
            call pic_gemm(current, rotation, updated)

            if (restricted) then
               call solve_ci(mol, updated, n_inactive, n_active, n_alpha, n_beta, &
                             alpha, beta, guess, have_guess, trial_ci, error, &
                             subspaces, min_electrons, max_electrons, flat_guess)
            else
               call solve_ci(mol, updated, n_inactive, n_active, n_alpha, n_beta, &
                             alpha, beta, guess, have_guess, trial_ci, error, &
                             n_roots=n_roots)
            end if
            if (error%has_error()) return
            trial_energy = trial_ci%energy
            if (averaged) trial_energy = weighted_sum(trial_ci%energies, weights)

            if (trial_energy < energy .or. predicted < ENERGY_RESOLUTION) then
               current = updated
               ci = trial_ci
               energy = trial_energy
               trust = min(MAX_ROTATION, trust*TRUST_GROWTH)
               accepted = .true.
               exit
            end if
            trust = 0.5_dp*trust
            if (trust < MIN_ROTATION) exit
         end do

         if (.not. accepted) then
            result%stalled = .true.
            call logger%warning("    no step downhill was found; stopping")
            exit
         end if
      end do

      result%energy = energy
      result%core_energy = ci%core_energy
      result%active_energy = ci%active_energy
      if (averaged) result%active_energy = energy - ci%core_energy
      result%gradient_norm = largest
      result%n_determinants = ci%n_determinants
      result%orbitals = current
      ! Whichever of the two the CI actually produced. A restricted space fills
      ! the flat one and leaves the rectangle unallocated, and copying an
      ! unallocated allocatable is undefined rather than empty.
      if (allocated(ci%ci_vector)) result%ci_vector = ci%ci_vector
      if (allocated(ci%ci_flat)) result%ci_flat = ci%ci_flat
      result%dm1 = dm1
      result%dm2 = dm2
      if (averaged) then
         allocate (result%energies(n_roots), result%spins(n_roots))
         result%energies = ci%energies(1:n_roots)
         result%ci_vectors = ci%vectors(:, :, 1:n_roots)
         do j = 1, n_roots
            result%spins(j) = spin_squared(n_active, n_alpha, n_beta, ci%vectors(:, :, j), &
                                           error)
            if (error%has_error()) return
         end do
      end if

      if (loud) then
         write (line, "(a,f22.12)") "    converged energy        ", result%energy
         call logger%info(trim(line))
         if (averaged) then
            do j = 1, n_roots
               write (line, "(a,i0,a,f20.12,a,f8.4)") "      state ", j, "  E = ", &
                  result%energies(j), "   <S^2> = ", result%spins(j)
               call logger%info(trim(line))
            end do
         end if
      end if

      ! Outside the `loud` block on purpose: a warning says the answer may be
      ! wrong, and that is not a thing a verbosity choice should be able to
      ! withhold. It was inside, so an MCSCF that never reached its gradient
      ! threshold said so only when someone had already asked for detail.
      if (.not. result%converged) then
         call logger%warning("    the orbital gradient did not reach the threshold")
      end if

      call alpha%destroy()
      call beta%destroy()
   end subroutine run_czt_casscf

   pure function weighted_sum(energies, weights) result(total)
      !! `sum_J weights(J) * energies(J)`, e.g. `E_SA` from every root's energy
      real(dp), intent(in) :: energies(:), weights(:)
      real(dp) :: total

      integer :: j

      total = 0.0_dp
      do j = 1, size(weights)
         total = total + weights(j)*energies(j)
      end do
   end function weighted_sum

   subroutine sa_density_matrices(vectors, weights, alpha, beta, dm1, dm2, error)
      !! `D_SA = sum_J w_J D_J`, `d_SA = sum_J w_J d_J`, one `active_space_rdms`
      !! call per state
      real(dp), intent(in) :: vectors(:, :, :)  !! (n_alpha_str, n_beta_str, n_states)
      real(dp), intent(in) :: weights(:)
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), allocatable, intent(out) :: dm1(:, :)
      real(dp), allocatable, intent(out) :: dm2(:, :, :, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dm1_j(:, :), dm2_j(:, :, :, :)
      integer :: j, n_states

      if (error%has_error()) return
      n_states = size(weights)
      do j = 1, n_states
         call active_space_rdms(vectors(:, :, j), alpha, beta, dm1_j, dm2_j, error)
         if (error%has_error()) return
         if (j == 1) then
            dm1 = weights(1)*dm1_j
            dm2 = weights(1)*dm2_j
         else
            dm1 = dm1 + weights(j)*dm1_j
            dm2 = dm2 + weights(j)*dm2_j
         end if
      end do
   end subroutine sa_density_matrices

   subroutine natural_orbitals(orbitals, n_inactive, n_active, dm1, natural, &
                               occupations, error)
      !! The orbitals that diagonalise the one-particle density, and their occupations
      !!
      !! An MCSCF wave function has no "occupied orbitals" in the sense a
      !! Hartree-Fock one does: the active orbitals carry fractional occupation
      !! and the optimised orbitals are not ordered by it at all. The natural
      !! orbitals are the closest thing -- the basis in which the density is
      !! diagonal -- and sorting them by occupation is what lets anything
      !! written for a reference determinant be pointed at a correlated one.
      !!
      !! The density is block diagonal (two on the inactive diagonal, the active
      !! block, zero on the virtual), so the inactive and virtual orbitals come
      !! through untouched. The whole matrix is diagonalised anyway, which means
      !! the Hartree-Fock case needs no separate path.
      real(dp), intent(in) :: orbitals(:, :)      !! (n_ao, n_mo)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1(:, :)           !! Active one-particle density
      real(dp), allocatable, intent(out) :: natural(:, :)
      real(dp), allocatable, intent(out) :: occupations(:)
         !! Descending, so the leading columns are the most occupied
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: density(:, :), values(:), ordered(:, :)
      integer :: n_ao, n_mo, i, info

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)

      allocate (density(n_mo, n_mo), values(n_mo))
      density = 0.0_dp
      do i = 1, n_inactive
         density(i, i) = 2.0_dp
      end do
      density(n_inactive + 1:n_inactive + n_active, &
              n_inactive + 1:n_inactive + n_active) = dm1

      call pic_syev(density, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the one-particle density could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
         return
      end if

      ! `pic_syev` returns ascending; occupations want the opposite.
      allocate (natural(n_ao, n_mo), occupations(n_mo), ordered(n_mo, n_mo))
      do i = 1, n_mo
         occupations(i) = values(n_mo - i + 1)
         ordered(:, i) = density(:, n_mo - i + 1)
      end do
      call pic_gemm(orbitals, ordered, natural)

      deallocate (density, values, ordered)
   end subroutine natural_orbitals

   subroutine solve_ci(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                       alpha, beta, guess, have_guess, ci, error, &
                       subspaces, min_electrons, max_electrons, flat_guess, n_roots)
      !! One CASCI, started from the previous vector when there is one
      !!
      !! After the first couple of macro-iterations the orbitals barely move, so
      !! the previous CI vector is very nearly the answer and the Davidson
      !! converges in a handful of products. A trust-region backtrack re-solves
      !! the CI at a rejected geometry, so a cheap re-solve is what makes
      !! rejecting a step affordable.
      type(czt_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), allocatable, intent(inout) :: guess(:, :, :)
      logical, intent(inout) :: have_guess
      type(casci_result_t), intent(out) :: ci
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: subspaces(:), min_electrons(:), max_electrons(:)
      real(dp), allocatable, intent(inout), optional :: flat_guess(:, :)
      integer, intent(in), optional :: n_roots
         !! More than one: that many singlet roots (state averaging), each kept
         !! as the next call's guess. Complete active space only.

      real(dp), parameter :: CI_TOLERANCE = 1.0e-11_dp
         !! Residual norm every CI solve here converges to, whichever branch

      ! A restricted space has no alpha-by-beta rectangle to keep a guess in, so
      ! it carries the flat vector instead.
      if (present(subspaces)) then
         if (have_guess .and. present(flat_guess)) then
            call run_czt_ormas_ci(mol, orbitals, n_inactive, n_active, n_alpha, &
                                  n_beta, subspaces, min_electrons, max_electrons, &
                                  ci, error, tolerance=CI_TOLERANCE, guess=flat_guess)
         else
            call run_czt_ormas_ci(mol, orbitals, n_inactive, n_active, n_alpha, &
                                  n_beta, subspaces, min_electrons, max_electrons, &
                                  ci, error, tolerance=CI_TOLERANCE)
         end if
         if (error%has_error()) return
         if (present(flat_guess)) then
            if (allocated(flat_guess)) deallocate (flat_guess)
            allocate (flat_guess(size(ci%ci_flat), 1))
            flat_guess(:, 1) = ci%ci_flat
            have_guess = .true.
         end if
         return
      end if

      if (present(n_roots)) then
         if (n_roots > 1) then
            if (have_guess) then
               call run_czt_casci(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                  ci, error, n_roots=n_roots, tolerance=CI_TOLERANCE, &
                                  guess=guess, symmetrize_singlet=.true.)
            else
               call run_czt_casci(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                  ci, error, n_roots=n_roots, tolerance=CI_TOLERANCE, &
                                  symmetrize_singlet=.true.)
            end if
            if (error%has_error()) return
            guess = ci%vectors
            have_guess = .true.
            return
         end if
      end if

      if (have_guess) then
         call run_czt_casci(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                            ci, error, tolerance=CI_TOLERANCE, guess=guess)
      else
         call run_czt_casci(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                            ci, error, tolerance=CI_TOLERANCE)
      end if
      if (error%has_error()) return

      if (allocated(guess)) deallocate (guess)
      allocate (guess(alpha%n_strings, beta%n_strings, 1))
      guess(:, :, 1) = ci%ci_vector
      have_guess = .true.
   end subroutine solve_ci

end module mqc_czt_mcscf
