!! The generalised Fock matrix, the orbital gradient and the orbital Hessian
module mqc_libcint_mcscf
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
   use pic_blas_interfaces, only: pic_gemm
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_libcint_integrals, only: libcint_molecule_t
   use mqc_libcint_direct, only: schwarz_bounds, build_fock_direct, direct_stats_t
   use mqc_libcint_mp2, only: transform_block
   use mqc_libcint_casci, only: run_libcint_casci, casci_result_t, &
                                run_libcint_ormas_ci
   use mqc_determinants, only: link_table_t, build_link_table
   use mqc_rdm, only: active_space_rdms
   use mqc_ormas_space, only: ormas_space_t, build_ormas_space
   use mqc_ormas_ci, only: ormas_density_matrices
   use pic_logger, only: logger => global_logger
   use pic_lapack_interfaces, only: pic_syev
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
   public :: run_libcint_casscf
   public :: casscf_result_t
   public :: natural_orbitals

   real(dp), parameter :: MIN_CURVATURE = 1.0e-3_dp
      !! Smallest curvature the Newton equations are allowed to divide by, in
      !! hartree per radian squared.
      !!
      !! Every eigenvalue is raised until the smallest reaches this, so a mode
      !! that is genuinely stiff keeps its own curvature and only the soft and
      !! the inverted ones are regularised. It has to sit below the softest real
      !! mode and above zero; larger damps directions that did not need it and
      !! costs the quadratic convergence the Hessian was built for.
   real(dp), parameter :: SADDLE_CURVATURE = -1.0e-5_dp
      !! How negative an eigenvalue has to be before the point is called a
      !! saddle rather than a minimum.
      !!
      !! Deliberately loose. A nearly redundant rotation -- an active orbital
      !! whose occupation has gone to two, say -- has a nearly zero eigenvalue
      !! whose sign is decided by how well the CI converged, and chasing those
      !! buys nothing. A real symmetry-breaking mode is orders of magnitude
      !! below this.
   real(dp), parameter :: MAX_ROTATION = 0.2_dp
      !! Largest rotation angle the trust radius is allowed to grow back to, in
      !! radians. Roughly 11 degrees.
   real(dp), parameter :: MIN_ROTATION = 1.0e-6_dp
      !! If backtracking has shrunk the trust radius this far and the energy
      !! still rises, the step direction is not a descent direction and halving
      !! it again will not help. Better to stop and say so.
   real(dp), parameter :: ENERGY_RESOLUTION = 1.0e-12_dp
      !! An energy change this small is not a change, in hartree.
      !!
      !! About fifteen units in the last place of a molecular energy, so it is
      !! the resolution of the arithmetic rather than a convergence criterion.
      !! Near the solution a Newton step is worth `g^2/H`, which drops below
      !! what the energy can report while the gradient is still shrinking, so a
      !! step whose *predicted* gain is below this is taken without being
      !! tested rather than rejected on the noise of the CI solve.
   real(dp), parameter :: TRUST_GROWTH = 1.3_dp
      !! How fast the trust radius recovers after a successful step. Slower than
      !! it shrinks: an over-long step costs a wasted CI solve, an over-short
      !! one only an iteration.

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
      ! TODO(mqc): calls `mol%eris_packed`, which recomputes the whole packed AO
      ! integral list rather than caching it. `orbital_hessian` does the same,
      ! so a CASSCF macro-iteration builds them twice.
      type(libcint_molecule_t), intent(in) :: mol
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

   subroutine rotation_matrix(kappa, rotation)
      !! `exp(kappa)` for antisymmetric `kappa`, by scaling and squaring
      !!
      !! Orthogonal to machine precision for any step size, so a large step is a
      !! bad step but never an invalid one and no reorthogonalisation is needed
      !! after it. Scaled by halving until the norm is below 1/2, because a
      !! plain Taylor series loses accuracy once the step is not small.
      real(dp), intent(in) :: kappa(:, :)
      real(dp), allocatable, intent(out) :: rotation(:, :)

      real(dp), allocatable :: scaled(:, :), term(:, :), next_term(:, :)
      real(dp) :: norm
      integer :: n, squarings, k, i
      integer, parameter :: TAYLOR_TERMS = 18

      n = size(kappa, 1)
      allocate (scaled(n, n), term(n, n), next_term(n, n), rotation(n, n))

      norm = maxval(abs(kappa))
      squarings = 0
      do while (norm > 0.5_dp)
         norm = 0.5_dp*norm
         squarings = squarings + 1
      end do
      scaled = kappa/real(2**squarings, dp)

      rotation = 0.0_dp
      term = 0.0_dp
      do i = 1, n
         rotation(i, i) = 1.0_dp
         term(i, i) = 1.0_dp
      end do
      do k = 1, TAYLOR_TERMS
         call pic_gemm(term, scaled, next_term)
         term = next_term/real(k, dp)
         rotation = rotation + term
         if (maxval(abs(term)) < 1.0e-18_dp) exit
      end do

      do k = 1, squarings
         call pic_gemm(rotation, rotation, next_term)
         rotation = next_term
      end do

      deallocate (scaled, term, next_term)
   end subroutine rotation_matrix

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
      real(dp), intent(in) :: a_block(:, :, :, :)   !! `(p q|r s)`, `q` and `s` occupied
      real(dp), intent(in) :: b_block(:, :, :, :)   !! `(p q|r s)`, `r` and `s` occupied
      integer, intent(in) :: n_occ
      real(dp), intent(in) :: density(:, :)         !! (n_mo, n_mo), symmetric
      real(dp), intent(out) :: potential(:, :)      !! (n_mo, n_occ)

      real(dp), allocatable :: weight(:, :), exchange(:, :)
      integer :: n_mo, q, r, s

      n_mo = size(density, 1)
      allocate (weight(n_mo, n_occ), exchange(n_mo, n_occ))

      ! The block on hand has only the second density index occupied. A pair
      ! with the *first* index occupied is the same integral read the other way
      ! round, so doubling the rows that are not occupied covers both.
      weight = density(:, 1:n_occ)
      weight(n_occ + 1:n_mo, :) = 2.0_dp*weight(n_occ + 1:n_mo, :)

      potential = 0.0_dp
      do s = 1, n_occ
         do r = 1, n_mo
            potential = potential + weight(r, s)*a_block(:, :, r, s)
         end do
      end do

      exchange = 0.0_dp
      do q = 1, n_occ
         do s = 1, n_occ
            do r = 1, n_mo
               exchange(:, q) = exchange(:, q) + density(r, s)*b_block(:, r, s, q)
            end do
         end do
         ! The exchange cannot be folded the same way: swapping the density
         ! indices moves them to different slots of the integral, so the half
         ! with a virtual second index comes from the other block.
         do s = n_occ + 1, n_mo
            do r = 1, n_occ
               exchange(:, q) = exchange(:, q) + density(r, s)*a_block(s, q, :, r)
            end do
         end do
      end do

      potential = potential - 0.5_dp*exchange
      deallocate (weight, exchange)
   end subroutine transformed_potential

   subroutine one_index_fock(a_block, b_block, fock, dm1, dm2, n_inactive, n_active, &
                             kappa, transformed)
      !! The generalised Fock matrix differentiated along one orbital rotation
      !!
      !!     (H kappa)_pq = 2 (Ft_qp - Ft_pq)
      !!
      !! with `Ft` built exactly as `generalized_fock` builds `F`, but with each
      !! orbital in turn replaced by `C kappa`. Differentiating an integral is
      !! the product rule over its coefficient matrices, so every term below is
      !! one of `generalized_fock`'s with one coefficient swapped.
      !!
      !! The density matrices are the ones the CI produced and are held fixed:
      !! this is the orbital-orbital block at fixed CI, not the full second
      !! derivative of the two-step energy. The difference is a Schur complement
      !! that only makes eigenvalues smaller, so a negative direction found here
      !! is a genuine one.
      real(dp), intent(in) :: a_block(:, :, :, :)   !! `(p q|r s)`, `q` and `s` occupied
      real(dp), intent(in) :: b_block(:, :, :, :)   !! `(p q|r s)`, `r` and `s` occupied
      type(mcscf_fock_t), intent(in) :: fock
      real(dp), intent(in) :: dm1(:, :), dm2(:, :, :, :)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: kappa(:, :)           !! (n_mo, n_mo), antisymmetric
      real(dp), intent(out) :: transformed(:, :)    !! (n_mo, n_mo)

      real(dp), allocatable :: d_inactive(:, :), d_active(:, :), occupations(:, :)
      real(dp), allocatable :: left(:, :), right(:, :)
      real(dp), allocatable :: v_inactive(:, :), v_active(:, :)
      real(dp), allocatable :: f_inactive(:, :), f_active(:, :), eri_gaaa(:, :, :, :)
      real(dp) :: accumulated
      integer :: n_mo, n_occ, i, t, u, v, w, n, ua, va, wa

      n_mo = size(kappa, 1)
      n_occ = n_inactive + n_active

      ! The inactive density is C_i C_i^T, so its derivative replaces one factor
      ! at a time and is non-zero only where exactly one index is inactive --
      ! which is also why a rotation between two inactive orbitals moves
      ! nothing.
      allocate (d_inactive(n_mo, n_mo), d_active(n_mo, n_mo), occupations(n_mo, n_mo))
      d_inactive = 0.0_dp
      d_inactive(:, 1:n_inactive) = 2.0_dp*kappa(:, 1:n_inactive)
      d_inactive(1:n_inactive, :) = d_inactive(1:n_inactive, :) &
                                    - 2.0_dp*kappa(1:n_inactive, :)

      occupations = 0.0_dp
      occupations(n_inactive + 1:n_occ, n_inactive + 1:n_occ) = dm1
      allocate (left(n_mo, n_mo), right(n_mo, n_mo))
      call pic_gemm(kappa, occupations, left)
      call pic_gemm(occupations, kappa, right)
      d_active = left - right

      allocate (v_inactive(n_mo, n_occ), v_active(n_mo, n_occ))
      call transformed_potential(a_block, b_block, n_occ, d_inactive, v_inactive)
      call transformed_potential(a_block, b_block, n_occ, d_active, v_active)

      ! Two contributions to each Fock matrix: the orbitals it is expressed in,
      ! which give a commutator, and the orbitals its density was built from,
      ! which give the potential above.
      allocate (f_inactive(n_mo, n_occ), f_active(n_mo, n_occ))
      call pic_gemm(fock%inactive, kappa, left)
      call pic_gemm(kappa, fock%inactive, right)
      f_inactive = left(:, 1:n_occ) - right(:, 1:n_occ) + v_inactive
      call pic_gemm(fock%active, kappa, left)
      call pic_gemm(kappa, fock%active, right)
      f_active = left(:, 1:n_occ) - right(:, 1:n_occ) + v_active

      ! (nu|vw) has four orbitals in it and so four terms, one per index.
      allocate (eri_gaaa(n_mo, n_active, n_active, n_active))
      do w = 1, n_active
         wa = n_inactive + w
         do v = 1, n_active
            va = n_inactive + v
            do u = 1, n_active
               ua = n_inactive + u
               do n = 1, n_mo
                  eri_gaaa(n, u, v, w) = &
                     dot_product(kappa(:, n), a_block(:, ua, va, wa)) &
                     + dot_product(kappa(:, ua), b_block(n, :, va, wa)) &
                     + dot_product(kappa(:, va), a_block(n, ua, :, wa)) &
                     + dot_product(kappa(:, wa), a_block(n, ua, :, va))
               end do
            end do
         end do
      end do

      ! The rows, assembled exactly as `generalized_fock` assembles them. The
      ! virtual rows stay zero because the density is, and the density is not
      ! what is being differentiated.
      transformed = 0.0_dp
      do i = 1, n_inactive
         do n = 1, n_mo
            transformed(i, n) = 2.0_dp*(f_inactive(n, i) + f_active(n, i))
         end do
      end do

      do t = 1, n_active
         do n = 1, n_mo
            accumulated = 0.0_dp
            do u = 1, n_active
               accumulated = accumulated + dm1(t, u)*f_inactive(n, n_inactive + u)
            end do
            do w = 1, n_active
               do v = 1, n_active
                  do u = 1, n_active
                     accumulated = accumulated + dm2(t, u, v, w)*eri_gaaa(n, u, v, w)
                  end do
               end do
            end do
            transformed(n_inactive + t, n) = accumulated
         end do
      end do

      deallocate (d_inactive, d_active, occupations, left, right)
      deallocate (v_inactive, v_active, f_inactive, f_active, eri_gaaa)
   end subroutine one_index_fock

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
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)      !! (n_ao, n_mo)
      integer, intent(in) :: n_inactive, n_active
      real(dp), intent(in) :: dm1(:, :), dm2(:, :, :, :)
      type(mcscf_fock_t), intent(in) :: fock
      integer, intent(in) :: rows(:), cols(:)     !! From `rotation_parameters`
      real(dp), allocatable, intent(out) :: hessian(:, :)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: eri_packed(:, :)
      real(dp), allocatable :: a_block(:, :, :, :), b_block(:, :, :, :)
      real(dp), allocatable :: kappa(:, :), transformed(:, :)
      integer :: n_mo, n_occ, n_param, k, l

      if (error%has_error()) return
      n_mo = size(orbitals, 2)
      n_occ = n_inactive + n_active
      n_param = size(rows)

      allocate (hessian(n_param, n_param))
      if (n_param == 0) return

      call mol%eris_packed(eri_packed)
      call transform_block(eri_packed, orbitals, orbitals(:, 1:n_occ), orbitals, &
                           orbitals(:, 1:n_occ), a_block)
      call transform_block(eri_packed, orbitals, orbitals, orbitals(:, 1:n_occ), &
                           orbitals(:, 1:n_occ), b_block)
      deallocate (eri_packed)

      allocate (kappa(n_mo, n_mo), transformed(n_mo, n_mo))
      do k = 1, n_param
         kappa = 0.0_dp
         kappa(rows(k), cols(k)) = 1.0_dp
         kappa(cols(k), rows(k)) = -1.0_dp
         call one_index_fock(a_block, b_block, fock, dm1, dm2, n_inactive, n_active, &
                             kappa, transformed)
         do l = 1, n_param
            hessian(l, k) = 2.0_dp*(transformed(cols(l), rows(l)) &
                                    - transformed(rows(l), cols(l)))
         end do
      end do

      ! Symmetric in exact arithmetic. Away from a stationary point the two
      ! triangles differ by rounding and by the antisymmetric term that
      ! distinguishes differentiating the gradient from differentiating the
      ! energy twice, which this average removes.
      hessian = 0.5_dp*(hessian + transpose(hessian))

      deallocate (a_block, b_block, kappa, transformed)
   end subroutine orbital_hessian

   subroutine newton_step(hessian, gradient, rows, cols, escape, kappa, lowest, &
                          predicted, error)
      !! The Newton step, level shifted, and pushed off a saddle when it is on one
      !!
      !! In the eigenbasis of the Hessian the step is one division per mode, and
      !! two things go wrong there.
      !!
      !! A mode with small or negative curvature would divide by nearly nothing,
      !! so every eigenvalue is raised by a single shift until the smallest
      !! reaches `MIN_CURVATURE`. One shift for all modes rather than a per-mode
      !! floor, because that is the exact solution of the trust-region
      !! subproblem and leaves the well-conditioned modes untouched.
      !!
      !! **A mode with negative curvature and no gradient on it is a saddle**,
      !! and it is where an optimiser built from the gradient alone stops and
      !! reports success. It happens whenever the starting orbitals carry a
      !! symmetry the solution does not. The division gives nothing to work with
      !! there, so such a mode is displaced by `escape` instead; either sign
      !! descends, and the caller backtracks if the step was too long.
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
         !! energy decrease. The caller uses it to tell a step that is too small
         !! to measure from one that is too long to take.
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: vectors(:, :), values(:), projected(:), amplitude(:)
      real(dp) :: shift, displacement
      integer :: n_mo, n_param, k, l, info

      lowest = 0.0_dp
      predicted = 0.0_dp
      if (error%has_error()) return
      n_mo = size(gradient, 1)
      n_param = size(rows)
      allocate (kappa(n_mo, n_mo))
      kappa = 0.0_dp
      if (n_param == 0) return

      allocate (vectors(n_param, n_param), values(n_param))
      allocate (projected(n_param), amplitude(n_param))
      vectors = hessian
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "the orbital Hessian could not be "// &
                        "diagonalized (info = "//to_char(info)//")")
         return
      end if
      lowest = values(1)

      shift = 0.0_dp
      if (values(1) < MIN_CURVATURE) shift = MIN_CURVATURE - values(1)

      do k = 1, n_param
         projected(k) = 0.0_dp
         do l = 1, n_param
            projected(k) = projected(k) + vectors(l, k)*gradient(rows(l), cols(l))
         end do
      end do

      do k = 1, n_param
         amplitude(k) = -projected(k)/(values(k) + shift)
         if (values(k) < SADDLE_CURVATURE .and. abs(amplitude(k)) < escape) then
            amplitude(k) = sign(escape, amplitude(k))
         end if
      end do

      do k = 1, n_param
         predicted = predicted - projected(k)*amplitude(k) &
                     - 0.5_dp*values(k)*amplitude(k)**2
      end do

      do l = 1, n_param
         displacement = dot_product(vectors(l, :), amplitude)
         kappa(rows(l), cols(l)) = displacement
         kappa(cols(l), rows(l)) = -displacement
      end do

      deallocate (vectors, values, projected, amplitude)
   end subroutine newton_step

   subroutine run_libcint_casscf(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                 result, error, max_iterations, gradient_tol, verbose, &
                                 subspaces, min_electrons, max_electrons)
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
      type(libcint_molecule_t), intent(in) :: mol
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

      type(ormas_space_t) :: space
      logical :: restricted
      type(casci_result_t) :: ci, trial_ci
      type(mcscf_fock_t) :: fock
      type(link_table_t) :: alpha, beta
      real(dp), allocatable :: current(:, :), updated(:, :)
      real(dp), allocatable :: dm1(:, :), dm2(:, :, :, :)
      real(dp), allocatable :: gradient(:, :), hessian(:, :), kappa(:, :), rotation(:, :)
      real(dp), allocatable :: step_matrix(:, :)
      real(dp), allocatable :: guess(:, :, :)
      real(dp), allocatable :: flat_guess(:, :)
      integer, allocatable :: rows(:), cols(:)
      character(len=160) :: line
      real(dp) :: tol, largest, previous, trust, scaling, lowest, predicted
      integer :: cycles, iteration, n_ao, n_mo, trial
      logical :: loud, have_guess, accepted, saddle

      integer, parameter :: MAX_BACKTRACKS = 12

      if (error%has_error()) return
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

      if (loud) then
         call logger%info("")
         if (restricted) then
            call logger%info("  occupation-restricted active space SCF")
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
                       guess, have_guess, ci, error)
      end if
      if (error%has_error()) return
      previous = ci%energy

      do iteration = 1, cycles
         result%iterations = iteration

         if (allocated(dm1)) deallocate (dm1)
         if (allocated(dm2)) deallocate (dm2)
         if (restricted) then
            call ormas_density_matrices(space, ci%ci_flat, dm1, dm2, error)
         else
            call active_space_rdms(ci%ci_vector, alpha, beta, dm1, dm2, error)
         end if
         if (error%has_error()) return

         call generalized_fock(mol, current, n_inactive, n_active, dm1, dm2, fock, error)
         if (error%has_error()) return
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
         call orbital_hessian(mol, current, n_inactive, n_active, dm1, dm2, fock, &
                              rows, cols, hessian, error)
         if (error%has_error()) return
         call newton_step(hessian, gradient, rows, cols, MAX_ROTATION, kappa, &
                          lowest, predicted, error)
         if (error%has_error()) return

         if (loud) then
            write (line, "(a,i4,f20.12,2es16.4,2es12.2)") "    ", iteration, &
               ci%energy, ci%energy - previous, largest, trust, lowest
            call logger%info(trim(line))
         end if
         previous = ci%energy

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
                             alpha, beta, guess, have_guess, trial_ci, error)
            end if
            if (error%has_error()) return

            if (trial_ci%energy < ci%energy .or. predicted < ENERGY_RESOLUTION) then
               current = updated
               ci = trial_ci
               trust = min(MAX_ROTATION, trust*TRUST_GROWTH)
               accepted = .true.
               exit
            end if
            trust = 0.5_dp*trust
            if (trust < MIN_ROTATION) exit
         end do

         if (.not. accepted) then
            result%stalled = .true.
            if (loud) call logger%warning("    no step downhill was found; stopping")
            exit
         end if
      end do

      result%energy = ci%energy
      result%core_energy = ci%core_energy
      result%active_energy = ci%active_energy
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

      if (loud) then
         write (line, "(a,f22.12)") "    converged energy        ", result%energy
         call logger%info(trim(line))
         if (.not. result%converged) then
            call logger%warning("    the orbital gradient did not reach the threshold")
         end if
      end if

      call alpha%destroy()
      call beta%destroy()
   end subroutine run_libcint_casscf

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
                       subspaces, min_electrons, max_electrons, flat_guess)
      !! One CASCI, started from the previous vector when there is one
      !!
      !! After the first couple of macro-iterations the orbitals barely move, so
      !! the previous CI vector is very nearly the answer and the Davidson
      !! converges in a handful of products. A trust-region backtrack re-solves
      !! the CI at a rejected geometry, so a cheap re-solve is what makes
      !! rejecting a step affordable.
      type(libcint_molecule_t), intent(in) :: mol
      real(dp), intent(in) :: orbitals(:, :)
      integer, intent(in) :: n_inactive, n_active, n_alpha, n_beta
      type(link_table_t), intent(in) :: alpha, beta
      real(dp), allocatable, intent(inout) :: guess(:, :, :)
      logical, intent(inout) :: have_guess
      type(casci_result_t), intent(out) :: ci
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: subspaces(:), min_electrons(:), max_electrons(:)
      real(dp), allocatable, intent(inout), optional :: flat_guess(:, :)

      ! A restricted space has no alpha-by-beta rectangle to keep a guess in, so
      ! it carries the flat vector instead.
      if (present(subspaces)) then
         if (have_guess .and. present(flat_guess)) then
            call run_libcint_ormas_ci(mol, orbitals, n_inactive, n_active, n_alpha, &
                                      n_beta, subspaces, min_electrons, max_electrons, &
                                      ci, error, tolerance=1.0e-11_dp, guess=flat_guess)
         else
            call run_libcint_ormas_ci(mol, orbitals, n_inactive, n_active, n_alpha, &
                                      n_beta, subspaces, min_electrons, max_electrons, &
                                      ci, error, tolerance=1.0e-11_dp)
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

      if (have_guess) then
         call run_libcint_casci(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                ci, error, tolerance=1.0e-11_dp, guess=guess)
      else
         call run_libcint_casci(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                ci, error, tolerance=1.0e-11_dp)
      end if
      if (error%has_error()) return

      if (allocated(guess)) deallocate (guess)
      allocate (guess(alpha%n_strings, beta%n_strings, 1))
      guess(:, :, 1) = ci%ci_vector
      have_guess = .true.
   end subroutine solve_ci

end module mqc_libcint_mcscf
