!! The generalised Fock matrix and the orbital gradient
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
   !! differences, not asserted.** Both orderings appear in the literature
   !! because both `C exp(kappa)` and `C exp(-kappa)` are used as the
   !! parametrisation, and getting it backwards gives an optimiser that climbs
   !! -- which looks like a convergence problem rather than a sign error, and is
   !! diagnosed as one for a long time. `test_mqc_mcscf.f90` differentiates the
   !! CASCI energy numerically and compares. Its rows are built differently
   !! depending on what the orbital is:
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
   use mqc_diis, only: diis_state_t
   implicit none
   private

   public :: mcscf_fock_t
   public :: generalized_fock
   public :: orbital_gradient
   public :: rotation_matrix
   public :: is_redundant
   public :: subspace_of
   public :: approximate_hessian
   public :: run_libcint_casscf
   public :: casscf_result_t
   public :: natural_orbitals

   real(dp), parameter :: HESSIAN_FLOOR = 0.05_dp
      !! Smallest denominator a step is allowed to divide by. Without it a pair
      !! of orbitals that happen to be close in energy takes an enormous step on
      !! a small gradient, which is how a first-order optimiser diverges on its
      !! first iteration rather than its fiftieth.
      !!
      !! **It is a genuine optimum, not a safety margin, and the direction is
      !! counter-intuitive.** The obvious reading of a slow case is that the
      !! floor is too high and the steps too short -- water crawls through a
      !! flat valley taking steps well inside the trust radius, which looks
      !! exactly like that. It is wrong. Lowering the floor makes it
      !! monotonically worse: 0.05, 0.02, 0.01, 0.005 and 0.002 cost 58, 87, 95,
      !! 145 and 192 iterations, because the longer steps overshoot, get
      !! rejected, and shrink the trust radius. Raising it is worse too -- 0.1,
      !! 0.2 and 0.4 cost 80, 115 and 206. Both directions were measured before
      !! this number was settled on.
   real(dp), parameter :: MAX_ROTATION = 0.2_dp
      !! Largest rotation angle the trust radius is allowed to grow back to, in
      !! radians. Roughly 11 degrees.
   real(dp), parameter :: MIN_ROTATION = 1.0e-6_dp
      !! If backtracking has shrunk the trust radius this far and the energy
      !! still rises, the step direction is not a descent direction and halving
      !! it again will not help. Better to stop and say so.
   integer, parameter :: DIIS_VECTORS = 12
      !! Orbital sets kept for extrapolation.
      !!
      !! Measured across 6, 10, 16 and 24 on water and nitrogen: nitrogen sits
      !! at 21 iterations regardless, and water moves between 56 and 75 with no
      !! trend -- which is the same spread it shows run to run at fixed
      !! settings, from the threaded Fock builds. So this is not a lever, and
      !! the value is a middle one rather than a tuned one.
   real(dp), parameter :: TRUST_GROWTH = 1.3_dp
      !! How fast the trust radius recovers after a successful step. Slower than
      !! it shrinks, which is the usual asymmetry: an over-long step costs a
      !! wasted CI solve, an over-short one only costs an iteration.

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
         !! rather than because it ran out of iterations. The distinction
         !! matters to whoever reads the failure: a stall is the approximate
         !! Hessian having run out of resolution, and more iterations do not
         !! help. It happens on flat surfaces -- water in cc-pVDZ with a
         !! CAS(4,4) stalls around a gradient of 1e-5 -- and the energy at that
         !! point is typically converged far past the gradient.
   end type casscf_result_t

   type :: mcscf_fock_t
      !! The generalised Fock matrix and the pieces it was built from
      real(dp), allocatable :: general(:, :)    !! (n_mo, n_mo), `F_mn`
      real(dp), allocatable :: inactive(:, :)   !! (n_mo, n_mo), `FI` in the MO basis
      real(dp), allocatable :: active(:, :)     !! (n_mo, n_mo), `FA` in the MO basis
      real(dp), allocatable :: occupation(:)
         !! (n_mo). Two for inactive, `D_tt` for active, zero for virtual. Used
         !! for the approximate Hessian, and worth having explicitly because it
         !! is what makes a rotation redundant when two orbitals share it.
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
      !! Getting this wrong is not loud. Treat a real parameter as redundant and
      !! the optimiser converges, reports a small gradient, and stops somewhere
      !! that is not a stationary point; treat a redundant one as real and the
      !! Hessian acquires a null direction. `subspaces` absent means a complete
      !! active space, which is the case the rest of this module was written
      !! for.
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

      ! The active mean field, which is the same J - K/2 built from the active
      ! density and with no one-electron part. Passing a zero core Hamiltonian
      ! is what leaves it out; there is no separate two-electron-only entry
      ! point and inventing one would be a second copy of the quartet loop.
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
      !! The result is orthogonal to machine precision for any step size, which
      !! is the reason for parametrising orbital changes this way: a large step
      !! is a bad step but never an invalid one, and no reorthogonalisation is
      !! needed after it.
      !!
      !! Scaling and squaring because a plain Taylor series loses accuracy once
      !! the step is not small, and orbital steps early in an optimisation are
      !! not small. Halving until the norm is below 1/2 costs a few matrix
      !! multiplies and makes the series converge in a handful of terms.
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

   subroutine approximate_hessian(fock, n_inactive, n_active, hessian, subspaces)
      !! A diagonal approximation to the orbital Hessian
      !!
      !! For a rotation between orbitals `p` and `q`,
      !!
      !!     H_pq ~ 2 |eps_p - eps_q| |n_p - n_q|
      !!
      !! with `eps` the diagonal of the combined inactive and active Fock and
      !! `n` the occupation. The two factors are the two ways a rotation can
      !! cost nothing: orbitals degenerate in energy, or equally occupied. Both
      !! are exactly the cases where the true Hessian is small or zero, so the
      !! approximation is crude in magnitude but right about where it vanishes,
      !! which is what a step needs it to be right about.
      !!
      !! Not the true Hessian and not intended to be. It scales the gradient
      !! into something the right order of magnitude; convergence is still
      !! first order and still needs acceleration. A second-order method would
      !! build the real thing, including the coupling between orbital rotations
      !! and the CI coefficients, and is a separate project.
      type(mcscf_fock_t), intent(in) :: fock
      integer, intent(in) :: n_inactive, n_active
      real(dp), allocatable, intent(out) :: hessian(:, :)
      integer, intent(in), optional :: subspaces(:)   !! As `orbital_gradient`

      real(dp), allocatable :: energies(:)
      integer :: n_mo, p, q

      n_mo = size(fock%general, 1)
      allocate (hessian(n_mo, n_mo), energies(n_mo))
      do p = 1, n_mo
         energies(p) = fock%inactive(p, p) + fock%active(p, p)
      end do

      hessian = HESSIAN_FLOOR
      do q = 1, n_mo
         do p = 1, n_mo
            if (is_redundant(p, q, n_inactive, n_active, subspaces)) cycle
            hessian(p, q) = max(HESSIAN_FLOOR, &
                                2.0_dp*abs(energies(p) - energies(q)) &
                                *abs(fock%occupation(p) - fock%occupation(q)))
         end do
      end do
      deallocate (energies)
   end subroutine approximate_hessian

   subroutine run_libcint_casscf(mol, orbitals, n_inactive, n_active, n_alpha, n_beta, &
                                 result, error, max_iterations, gradient_tol, verbose, &
                                 subspaces, min_electrons, max_electrons)
      !! Two-step CASSCF: solve the CI, move the orbitals, repeat
      !!
      !! Each macro-iteration solves the CI problem exactly at the current
      !! orbitals, builds the density matrices from it, and takes one orbital
      !! step from the gradient. The CI is re-solved from the previous vector,
      !! which after the first few iterations is nearly the answer already.
      !!
      !! **Two-step, so convergence is first order.** The orbital step ignores
      !! the coupling between orbital rotations and the CI coefficients -- it
      !! assumes the CI is re-optimised afterwards, which it is, but not that
      !! the two respond to each other. A second-order method treats them
      !! together and converges quadratically. This is perhaps a fifth of the
      !! work and gets to a validated energy, after which a better optimiser is
      !! a change of step rule against a known answer rather than a new method
      !! with two places to be wrong.
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
      real(dp), allocatable :: step_matrix(:, :), overlap(:, :), flat(:)
      real(dp), allocatable :: candidate(:, :)
      real(dp), allocatable :: guess(:, :, :)
      real(dp), allocatable :: flat_guess(:, :)
      type(diis_state_t) :: diis
      logical :: extrapolated
      character(len=160) :: line
      real(dp) :: tol, largest, previous, trust, scaling
      integer :: cycles, iteration, n_ao, n_mo, trial, diis_failures
      logical :: loud, have_guess, accepted

      integer, parameter :: MAX_BACKTRACKS = 12

      if (error%has_error()) return
      cycles = 50
      if (present(max_iterations)) cycles = max_iterations
      tol = 1.0e-6_dp
      if (present(gradient_tol)) tol = gradient_tol
      loud = .false.
      if (present(verbose)) loud = verbose

      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)
      allocate (current(n_ao, n_mo), updated(n_ao, n_mo), step_matrix(n_mo, n_mo))
      allocate (candidate(n_ao, n_mo))
      diis_failures = 0
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

      call mol%overlap(overlap)
      call diis%init(DIIS_VECTORS, n_ao*n_mo, n_mo*n_mo)
      allocate (flat(n_ao*n_mo))

      if (loud) then
         call logger%info("")
         if (restricted) then
            call logger%info("  occupation-restricted active space SCF")
         else
            call logger%info("  complete active space SCF")
         end if
         call logger%info("    iter            energy          change      max gradient"// &
                          "       trust")
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

         if (loud) then
            write (line, "(a,i4,f20.12,2es16.4,es12.2)") "    ", iteration, ci%energy, &
               ci%energy - previous, largest, trust
            call logger%info(trim(line))
         end if
         previous = ci%energy

         if (largest < tol) then
            result%converged = .true.
            exit
         end if

         if (allocated(hessian)) deallocate (hessian)
         if (restricted) then
            call approximate_hessian(fock, n_inactive, n_active, hessian, subspaces)
         else
            call approximate_hessian(fock, n_inactive, n_active, hessian)
         end if
         if (allocated(kappa)) deallocate (kappa)
         allocate (kappa(n_mo, n_mo))
         kappa = -gradient/hessian

         ! **The trust region is what makes this converge rather than
         ! oscillate.** The approximate Hessian is crude enough that the step it
         ! suggests is regularly too long near the solution, and a first-order
         ! method with a fixed step length does not settle -- it bounces across
         ! the minimum indefinitely, which is exactly what this did before the
         ! backtracking was added. Rejecting any step that raises the energy and
         ! halving the radius costs a CI solve when it happens and removes the
         ! failure mode entirely.
         accepted = .false.

         ! ---- extrapolate, once, at full step length -----------------------
         !
         ! **A failed extrapolation must not shrink the trust radius.** The
         ! trust radius governs the *gradient* step, and DIIS failing says
         ! nothing about whether that step was too long -- it says the subspace
         ! has gone stale. Conflating the two was measured: the trust radius
         ! collapsed to 3e-3 against a gradient of 2e-2, every step became
         ! useless, and water stalled indefinitely at -76.042 having converged
         ! cleanly without DIIS at all. So the extrapolation gets one attempt of
         ! its own, and only the gradient step below can shrink the radius.
         scaling = maxval(abs(kappa))
         if (scaling > trust) then
            step_matrix = kappa*(trust/scaling)
         else
            step_matrix = kappa
         end if
         if (allocated(rotation)) deallocate (rotation)
         call rotation_matrix(step_matrix, rotation)
         call pic_gemm(current, rotation, updated)

         call diis%push(reshape(updated, [n_ao*n_mo]), reshape(gradient, [n_mo*n_mo]))
         if (iteration > 1) then
            flat = reshape(updated, [n_ao*n_mo])
            call diis%extrapolate(flat, extrapolated)
            if (extrapolated) then
               candidate = reshape(flat, [n_ao, n_mo])
               call symmetric_orthonormalize(overlap, candidate, error)
               if (error%has_error()) return
               if (restricted) then
                  call solve_ci(mol, candidate, n_inactive, n_active, n_alpha, n_beta, &
                                alpha, beta, guess, have_guess, trial_ci, error, &
                                subspaces, min_electrons, max_electrons, flat_guess)
               else
                  call solve_ci(mol, candidate, n_inactive, n_active, n_alpha, n_beta, &
                                alpha, beta, guess, have_guess, trial_ci, error)
               end if
               if (error%has_error()) return
               if (trial_ci%energy < ci%energy) then
                  current = candidate
                  ci = trial_ci
                  accepted = .true.
                  diis_failures = 0
               else
                  ! A stale subspace keeps producing the same bad direction, so
                  ! it is cleared rather than left to fail every iteration.
                  diis_failures = diis_failures + 1
                  if (diis_failures >= 2) then
                     call diis%init(DIIS_VECTORS, n_ao*n_mo, n_mo*n_mo)
                     diis_failures = 0
                  end if
               end if
            end if
         end if

         ! ---- otherwise the gradient step, backtracking until it descends ---
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

            if (trial_ci%energy < ci%energy) then
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
      ! unallocated allocatable is undefined rather than empty -- it survived
      ! one set of compiler flags and segfaulted under another.
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

   subroutine symmetric_orthonormalize(overlap, orbitals, error)
      !! Restore orthonormality to orbitals that a linear combination broke
      !!
      !! `C <- C (C^T S C)^{-1/2}`, which is the change that moves every orbital
      !! as little as possible -- the right choice here because the input is
      !! already nearly orthonormal and the aim is to repair it, not to impose
      !! some other structure on it.
      !!
      !! Needed because DIIS extrapolates a *linear* combination of orbital sets
      !! and orthonormality is not linear. This is the price of extrapolating
      !! the orbitals directly instead of the rotation parameters that produced
      !! them; the alternative needs a matrix logarithm to recover those
      !! parameters from a product of exponentials, which is more machinery for
      !! the same answer.
      real(dp), intent(in) :: overlap(:, :)      !! AO overlap
      real(dp), intent(inout) :: orbitals(:, :)  !! (n_ao, n_mo)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: metric(:, :), vectors(:, :), values(:)
      real(dp), allocatable :: scaled(:, :), inverse_root(:, :), work(:, :)
      integer :: n_ao, n_mo, i, info

      if (error%has_error()) return
      n_ao = size(orbitals, 1)
      n_mo = size(orbitals, 2)

      allocate (work(n_ao, n_mo), metric(n_mo, n_mo))
      call pic_gemm(overlap, orbitals, work)
      call pic_gemm(orbitals, work, metric, transa="T")

      allocate (vectors(n_mo, n_mo), values(n_mo), scaled(n_mo, n_mo))
      allocate (inverse_root(n_mo, n_mo))
      vectors = metric
      call pic_syev(vectors, values, jobz="V", uplo="U", info=info)
      if (info /= 0 .or. values(1) <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "the extrapolated orbitals are linearly "// &
                        "dependent and cannot be orthonormalized. The DIIS subspace "// &
                        "has produced a combination that is not a set of orbitals.")
         return
      end if
      do i = 1, n_mo
         scaled(:, i) = vectors(:, i)/sqrt(values(i))
      end do
      call pic_gemm(scaled, vectors, inverse_root, transb="T")
      call pic_gemm(orbitals, inverse_root, work)
      orbitals = work

      deallocate (work, metric, vectors, values, scaled, inverse_root)
   end subroutine symmetric_orthonormalize

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
      !! through untouched and only the active block is rotated. Diagonalising
      !! the whole thing anyway costs nothing at this size and means the
      !! Hartree-Fock case needs no separate path: there the density is already
      !! diagonal and the "rotation" is the identity.
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
      !! converges in a handful of products instead of dozens. That matters
      !! here more than anywhere: a trust-region backtrack solves the CI again
      !! at a rejected geometry, so a cheap re-solve is what makes rejecting a
      !! step affordable.
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
      ! it carries the flat vector instead. Same idea and the same payoff: after
      ! the first few macro-iterations the previous answer is nearly this one,
      ! and a trust-region backtrack re-solves the CI at a rejected geometry.
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
