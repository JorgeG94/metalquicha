module mqc_scf_gdm
   !! Numerical core of geometric direct minimisation (GDM) for the SCF.
   !!
   !! GDM (Van Voorhis and Head-Gordon, Mol. Phys. 100, 1713 (2002)) minimises
   !! the SCF energy directly over the occupied-virtual orbital rotations:
   !! C <- C exp(kappa), with kappa chosen along a quasi-Newton direction and
   !! a step accepted only if the energy went down. This module holds the
   !! pieces of that which involve no SCF, basis or device code, so they can be
   !! tested against dense linear algebra on their own:
   !!
   !! - `scf_gdm_history_t`: a limited-memory BFGS history (Nocedal, Math.
   !!   Comp. 35, 773 (1980)) whose two-loop recursion applies the inverse
   !!   Hessian to the gradient without ever forming it;
   !! - `scf_gdm_preconditioner_build`: the diagonal initial inverse Hessian
   !!   from the orbital energies (the preconditioner);
   !! - `scf_gdm_step_judge`: the accept-only-on-energy-decrease rule and the
   !!   step scale it backs off with;
   !! - `scf_gdm_step_limit`: the factor that caps the largest rotation;
   !! - `scf_gdm_is_descent`: the check that a direction goes downhill.
   !!
   !! MISUSE. A wrong vector length, an uninitialised history or a bad
   !! parameter is a programming error, reported through the `err` argument
   !! (`mqc_error`, code `ERROR_VALIDATION`) and never confused with the
   !! ordinary outcomes (a pair that fails the curvature test comes back as
   !! `accepted = .false.` with `err` untouched). Each such procedure leaves
   !! its outputs untouched on misuse, and the caller clears `err` first.
   !!
   !! VECTOR LAYOUT. A rotation vector of one spin block has length nocc*nvir
   !! and element (a, i) at position (i-1)*nvir + a: the virtual index `a`
   !! runs fastest, the occupied index `i` slowest. `scf_gdm_preconditioner_build`
   !! produces the preconditioner in that order and the SCF step must pack its
   !! gradient and rotation amplitudes the same way. An unrestricted run
   !! concatenates the alpha block and then the beta block into one vector.
   !!
   !! Plain arrays throughout; the vectors are at most a few million doubles.
   !! The same module, name for name, as terco's `trc_scf_gdm`, so the two
   !! backends take the same step from the same history.
   use pic_types, only: dp
   use pic_blas_interfaces, only: pic_dot
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: scf_gdm_history_t
   public :: scf_gdm_preconditioner_build, scf_gdm_step_judge, scf_gdm_step_limit, scf_gdm_is_descent
   public :: scf_gdm_curvature_tol, scf_gdm_scale_min, scf_gdm_descent_margin

   real(dp), parameter :: scf_gdm_curvature_tol = 1.0e-4_dp
      !! A pair (s, y) is stored only if the cosine s.y/(|s||y|) exceeds this.
      !! BFGS keeps the inverse Hessian positive definite exactly when s.y > 0,
      !! but a pair with a tiny positive cosine is noise: with H0 = I and
      !! |s| = |y| the direction grows like 1/cosine^2 (measured |d|/|g| = 5.7e4
      !! at cos = 2e-3 and 1.6e18 at 2e-10), more than the step scale
      !! back-off can undo; `scf_gdm_step_limit` caps what remains. A legitimate pair on a quadratic of condition
      !! number kappa has cosine >= 2 sqrt(kappa)/(kappa + 1), i.e. 6e-4 at
      !! kappa = 1e7, so 1e-4 loses no real curvature information.
   real(dp), parameter :: scf_gdm_scale_min = 1.0_dp/1024.0_dp
      !! The step scale below which `scf_gdm_step_judge` reports breakdown:
      !! a scale of 1 halved eleven times in a row (2^-11 < 2^-10) without an
      !! energy decrease means the direction is not a usable one and the
      !! caller must restart from the last accepted point or stop.
   real(dp), parameter :: scf_gdm_descent_margin = 1.0e-12_dp
      !! `scf_gdm_is_descent` needs g.d < -margin |g||d|: a direction that is
      !! orthogonal to the gradient to rounding is not a descent direction.

   type :: scf_gdm_history_t
      !! Circular buffer of the last `m` accepted pairs s_k = x_{k+1} - x_k
      !! (the step) and y_k = g_{k+1} - g_k (the gradient change). Declared
      !! with plain `type ::` and filled by `init`, never by a structure
      !! constructor (allocatable components, their defaults do the rest).
      integer :: n = 0
         !! Length of every vector.
      integer :: m = 0
         !! Capacity: the most pairs kept; older ones are overwritten.
      integer :: count = 0
         !! Pairs currently stored (at most m).
      integer :: head = 1
         !! Slot the next accepted pair is written to.
      real(dp), allocatable :: s(:, :)
         !! Step vectors, shape (n, m).
      real(dp), allocatable :: y(:, :)
         !! Gradient-change vectors, shape (n, m).
      real(dp), allocatable :: rho(:)
         !! 1/(s_k.y_k) per slot, kept so the two loops do not redo the dot.
   contains
      procedure :: init => scf_gdm_history_init
      procedure :: reset => scf_gdm_history_reset
      procedure :: release => scf_gdm_history_release
      procedure :: push => scf_gdm_history_push
      procedure :: direction => scf_gdm_history_direction
      procedure :: npairs => scf_gdm_history_npairs
      procedure, private :: slot => scf_gdm_history_slot
   end type scf_gdm_history_t

contains

   subroutine scf_gdm_history_init(self, n, m, err)
      !! Allocate room for `m` pairs of length `n` and start empty. Calling it
      !! on an initialised history releases the old storage first. The pair
      !! storage is deliberately left untouched after allocation (2*n*m
      !! doubles, ~1 GB at n = 3e6, m = 20): only stored slots are ever read.
      class(scf_gdm_history_t), intent(inout) :: self
      integer, intent(in) :: n
         !! Vector length, at least 1.
      integer, intent(in) :: m
         !! History capacity, at least 1.
      type(error_t), intent(inout) :: err
         !! Set if n < 1 or m < 1; the history is then left released.
      call self%release()
      if (n < 1 .or. m < 1) then
         call err%set(ERROR_VALIDATION, "scf_gdm: history needs n >= 1 and m >= 1")
         return
      end if
      self%n = n
      self%m = m
      allocate (self%s(n, m), self%y(n, m), self%rho(m))
      self%rho = 0.0_dp
      self%count = 0
      self%head = 1
   end subroutine scf_gdm_history_init

   subroutine scf_gdm_history_reset(self)
      !! Forget every pair but keep the allocation. The SCF step does this when
      !! the quasi-Newton direction fails the descent test, so the next step
      !! is plain preconditioned steepest descent.
      class(scf_gdm_history_t), intent(inout) :: self
      self%count = 0
      self%head = 1
   end subroutine scf_gdm_history_reset

   subroutine scf_gdm_history_release(self)
      !! Free the storage; the history is then uninitialised (n = m = 0).
      class(scf_gdm_history_t), intent(inout) :: self
      if (allocated(self%s)) deallocate (self%s)
      if (allocated(self%y)) deallocate (self%y)
      if (allocated(self%rho)) deallocate (self%rho)
      self%n = 0
      self%m = 0
      self%count = 0
      self%head = 1
   end subroutine scf_gdm_history_release

   pure integer function scf_gdm_history_npairs(self) result(npairs)
      !! Number of pairs currently stored.
      class(scf_gdm_history_t), intent(in) :: self
      npairs = self%count
   end function scf_gdm_history_npairs

   pure integer function scf_gdm_history_slot(self, age) result(slot)
      !! Buffer slot of the pair of the given age, 1 = oldest, count = newest.
      class(scf_gdm_history_t), intent(in) :: self
      integer, intent(in) :: age
      slot = modulo(self%head - 1 - self%count + age - 1, self%m) + 1
   end function scf_gdm_history_slot

   subroutine scf_gdm_history_push(self, s, y, accepted, err)
      !! Store the pair (s, y) unless it would spoil the BFGS update.
      !!
      !! The BFGS inverse Hessian stays positive definite if and only if every
      !! stored pair has s.y > 0. A rotation step on a non-convex energy can
      !! produce s.y <= 0, and a failed or noisy step can produce NaN; such a
      !! pair is dropped (`accepted` = .false.) and the history is untouched,
      !! so the direction afterwards is still a descent direction. The test is
      !! written so a NaN anywhere fails it. An uninitialised history or a
      !! vector of the wrong length is a programming error, reported in `err`
      !! (and `accepted` is false), not a rejection.
      class(scf_gdm_history_t), intent(inout) :: self
      real(dp), intent(in) :: s(:)
         !! Step taken, length n.
      real(dp), intent(in) :: y(:)
         !! Change of the gradient over that step, length n.
      logical, intent(out) :: accepted
         !! True if the pair was stored.
      type(error_t), intent(inout) :: err
         !! Set on an uninitialised history or a length mismatch.
      real(dp) :: sy, ss, yy
      integer :: slot

      accepted = .false.
      if (self%m < 1 .or. size(s) /= self%n .or. size(y) /= self%n) then
         call err%set(ERROR_VALIDATION, "scf_gdm: push on an uninitialised history or with a wrong vector length")
         return
      end if
      sy = pic_dot(s, y)
      ss = pic_dot(s, s)
      yy = pic_dot(y, y)
      ! Inf or NaN in an element reaches all three dots as Inf or NaN, and
      ! every comparison below is false for NaN, so no separate scan is needed.
      if (.not. (ss > 0.0_dp .and. yy > 0.0_dp .and. ss <= huge(1.0_dp) &
                 .and. yy <= huge(1.0_dp) .and. abs(sy) <= huge(1.0_dp))) return
      if (.not. (sy > scf_gdm_curvature_tol*sqrt(ss)*sqrt(yy))) return
      ! sqrt(ss)*sqrt(yy) can still overflow or underflow for huge or tiny
      ! vectors; the comparison then rejects, which is the safe side.

      slot = self%head
      self%s(:, slot) = s
      self%y(:, slot) = y
      self%rho(slot) = 1.0_dp/sy
      self%head = modulo(slot, self%m) + 1
      self%count = min(self%count + 1, self%m)
      accepted = .true.
   end subroutine scf_gdm_history_push

   subroutine scf_gdm_history_direction(self, g, h0_inv, d, err)
      !! d = -H g, with H the L-BFGS inverse Hessian built from the stored
      !! pairs and the DIAGONAL initial inverse Hessian `h0_inv`.
      !!
      !! This is the two-loop recursion of Nocedal (1980). The first loop walks
      !! the pairs newest to oldest, stripping from q = g the part of the
      !! gradient each pair explains; the diagonal then plays the role of the
      !! Hessian that is not covered by any pair; the second loop walks oldest
      !! to newest and puts the pair corrections back. Cost is O(m n), no
      !! matrix is formed. The diagonal is the preconditioner (about the
      !! inverse of the orbital-energy gap) and already carries the scale of
      !! the problem, so there is deliberately no Shanno gamma = s.y/y.y
      !! rescaling on top of it. With no pairs stored it returns the
      !! preconditioned steepest-descent direction -h0_inv*g.
      class(scf_gdm_history_t), intent(in) :: self
      real(dp), intent(in) :: g(:)
         !! Gradient, length n.
      real(dp), intent(in) :: h0_inv(:)
         !! Diagonal of the initial inverse Hessian, length n, positive.
      real(dp), intent(out) :: d(:)
         !! Search direction, length n; untouched if `err` is set.
      type(error_t), intent(inout) :: err
         !! Set on an uninitialised history or any length that is not n.
      real(dp) :: alpha(self%count)
      real(dp) :: beta
      integer :: age, slot

      if (self%m < 1 .or. size(g) /= self%n .or. size(h0_inv) /= self%n .or. size(d) /= self%n) then
         call err%set(ERROR_VALIDATION, "scf_gdm: direction on an uninitialised history or with a wrong vector length")
         return
      end if
      d = g
      do age = self%count, 1, -1
         slot = self%slot(age)
         alpha(age) = self%rho(slot)*pic_dot(self%s(:, slot), d)
         d = d - alpha(age)*self%y(:, slot)
      end do
      d = h0_inv*d
      do age = 1, self%count
         slot = self%slot(age)
         beta = self%rho(slot)*pic_dot(self%y(:, slot), d)
         d = d + (alpha(age) - beta)*self%s(:, slot)
      end do
      d = -d
   end subroutine scf_gdm_history_direction

   pure subroutine scf_gdm_preconditioner_build(e_occ, e_vir, hessian_factor, hessian_floor, h0_inv, err)
      !! Diagonal initial inverse Hessian for one spin block.
      !!
      !! To leading order the second derivative of the energy with respect to
      !! the rotation kappa_ai is the Hessian estimate
      !! hessian_factor*(e_a - e_i), where `hessian_factor` is the factor
      !! between dE/dkappa_ai and F_ai (4 for a restricted closed shell, 2 per
      !! spin when unrestricted). Its inverse is the natural scale of a step
      !! in that coordinate. The ESTIMATE hessian_factor*gap (not the bare
      !! orbital-energy gap) is floored at `hessian_floor`, so a tiny or
      !! negative gap (a near-degenerate or level-crossing orbital pair)
      !! cannot make a step enormous or flip its sign; the floor is in units
      !! of the Hessian estimate. A non-positive factor or floor, or an
      !! `h0_inv` whose length is not nocc*nvir, is a programming error
      !! reported in `err` (`h0_inv` untouched).
      !!
      !! Index convention: h0_inv((i-1)*nvir + a) belongs to the rotation of
      !! occupied i into virtual a; `a` is the fast index, `i` the slow one.
      real(dp), intent(in) :: e_occ(:)
         !! Occupied diagonal Fock elements, length nocc.
      real(dp), intent(in) :: e_vir(:)
         !! Virtual diagonal Fock elements, length nvir.
      real(dp), intent(in) :: hessian_factor
         !! 4 (restricted closed shell) or 2 (per spin, unrestricted), positive.
      real(dp), intent(in) :: hessian_floor
         !! Smallest Hessian estimate used (units of the estimate), positive.
      real(dp), intent(out) :: h0_inv(:)
         !! Preconditioner, length nocc*nvir.
      type(error_t), intent(inout) :: err
         !! Set on a bad factor, floor or output length.
      integer :: i, a, nvir

      nvir = size(e_vir)
      if (.not. (hessian_factor > 0.0_dp .and. hessian_floor > 0.0_dp) .or. size(h0_inv) /= size(e_occ)*nvir) then
         call err%set(ERROR_VALIDATION, "scf_gdm: preconditioner needs a positive factor and floor and nocc*nvir outputs")
         return
      end if
      do i = 1, size(e_occ)
         do a = 1, nvir
            h0_inv((i - 1)*nvir + a) = 1.0_dp/max(hessian_factor*(e_vir(a) - e_occ(i)), hessian_floor)
         end do
      end do
   end subroutine scf_gdm_preconditioner_build

   pure subroutine scf_gdm_step_judge(energy_last, energy_trial, slack, step_scale, accepted, broken)
      !! Accept-only-on-decrease rule and the step-scale back-off.
      !!
      !! The trial is accepted when it is finite and no higher than the last
      !! accepted energy plus `slack` (the noise floor of the energy, so a
      !! converged step that changes E by rounding is not thrown away). On
      !! acceptance the scale doubles back toward 1; on rejection it halves,
      !! and once it would fall below `scf_gdm_scale_min` (on the 11th
      !! consecutive rejection from 1) the routine sets `broken` so the caller
      !! stops retrying instead of shrinking a step forever. On breakdown
      !! `step_scale` is returned as 1 so a restart does not break down again
      !! after one rejection: the caller must go back to the last accepted
      !! point, reset the history and build a fresh direction. A NaN or Inf
      !! trial energy is always a rejection.
      real(dp), intent(in) :: energy_last
         !! Last accepted energy.
      real(dp), intent(in) :: energy_trial
         !! Energy at the trial orbitals.
      real(dp), intent(in) :: slack
         !! Non-negative tolerance on an increase.
      real(dp), intent(inout) :: step_scale
         !! Current step scale in (0, 1]; updated for the next trial, 1 on `broken`.
      logical, intent(out) :: accepted
         !! True if the trial is kept.
      logical, intent(out) :: broken
         !! True if the scale fell below the minimum: give up on this step.
      broken = .false.
      accepted = energy_trial <= energy_last + slack
      if (accepted .and. abs(energy_trial) > huge(1.0_dp)) accepted = .false.
      if (accepted) then
         step_scale = min(1.0_dp, 2.0_dp*step_scale)
      else
         step_scale = 0.5_dp*step_scale
         broken = step_scale < scf_gdm_scale_min
         if (broken) step_scale = 1.0_dp
      end if
   end subroutine scf_gdm_step_judge

   pure subroutine scf_gdm_step_limit(d, max_rotation, factor, valid)
      !! The factor in (0, 1] by which a direction must be scaled so that its
      !! largest element is at most `max_rotation` (an orbital rotation angle
      !! beyond about 1 radian leaves the region where the energy is
      !! predictable). 1 when it already fits, for a zero-length or all-zero
      !! `d`. A non-finite element in `d`, or a `max_rotation` that is not
      !! positive and finite, sets `valid` = .false. and `factor` = 1: the
      !! caller must reject the direction, not scale it.
      real(dp), intent(in) :: d(:)
         !! Direction.
      real(dp), intent(in) :: max_rotation
         !! Largest allowed |d(k)|, positive.
      real(dp), intent(out) :: factor
         !! Scale for `d`.
      logical, intent(out) :: valid
         !! False on a non-finite `d` or an unusable `max_rotation`.
      real(dp) :: largest

      factor = 1.0_dp
      valid = max_rotation > 0.0_dp .and. max_rotation <= huge(1.0_dp)
      if (size(d) == 0) return
      largest = maxval(abs(d))
      ! A NaN makes the comparison false, so it is caught here too.
      if (.not. (all(abs(d) <= huge(1.0_dp)))) valid = .false.
      if (.not. valid .or. largest <= max_rotation) return
      factor = max_rotation/largest
      ! Division can round up by one ulp; the contract is a hard cap.
      if (factor*largest > max_rotation) factor = nearest(factor, -1.0_dp)
   end subroutine scf_gdm_step_limit

   logical function scf_gdm_is_descent(g, d) result(is_descent)
      !! True if g.d < -scf_gdm_descent_margin |g||d|, i.e. moving along d
      !! lowers the energy to first order by more than rounding. The SCF step
      !! falls back to preconditioned steepest descent and resets the history
      !! when this fails. NaN gives false.
      real(dp), intent(in) :: g(:)
         !! Gradient.
      real(dp), intent(in) :: d(:)
         !! Proposed direction.
      is_descent = pic_dot(g, d) < -scf_gdm_descent_margin*sqrt(pic_dot(g, g))*sqrt(pic_dot(d, d))
   end function scf_gdm_is_descent

end module mqc_scf_gdm
