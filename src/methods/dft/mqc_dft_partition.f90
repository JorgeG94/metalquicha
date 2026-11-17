module mqc_dft_partition
   !! Becke partition of space into atomic cells
   !!
   !! A molecular integral is done as a sum of atom-centred quadratures, which
   !! requires a partition of unity: a set of weights w_A(r) summing to 1 at
   !! every point, each concentrated near its atom. Becke's construction
   !! [JCP 88, 2547 (1988)] builds them from the confocal elliptical coordinate
   !! of each atom pair,
   !!
   !!    mu_AB = (|r - R_A| - |r - R_B|) / R_AB
   !!
   !! smoothed by a cutoff s(mu) that runs from 1 to 0 as mu goes -1 to 1. The
   !! unnormalised cell function is the product over partners,
   !!
   !!    P_A(r) = prod_{B /= A} s(nu_AB),   w_A = P_A / sum_C P_C
   !!
   !! Two cutoffs are offered. Becke's is three iterations of p(x) = x(3-x^2)/2,
   !! smooth and non-zero everywhere, so every atom contributes at every point.
   !! Stratmann's [CPL 257, 213 (1996)] is a degree-5 polynomial reaching
   !! exactly +-1 at |mu| >= a, so a distant pair contributes exactly 0 or 1 and
   !! could be skipped.
   !!
   !! `nu` is `mu` after an adjustment for the two atoms having different sizes,
   !! so that the cell boundary sits between them rather than halfway. Becke's
   !! appendix gives the shift in terms of the radius ratio chi:
   !!
   !!    u = (chi - 1)/(chi + 1),   a = u/(u^2 - 1),  clamped to |a| <= 1/2
   !!
   !! which reduces to a = (1/chi - chi)/4. Treutler's variant uses the square
   !! roots of the radii instead. Both are available; the choice changes the
   !! weights slightly and must match whatever a reference is being compared to.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_dft_radial, only: bragg_radius
   implicit none
   private

   public :: PARTITION_BECKE, PARTITION_STRATMANN
   public :: ADJUST_NONE, ADJUST_BECKE, ADJUST_TREUTLER
   public :: becke_partition_weights
   public :: becke_partition_derivatives
   public :: becke_cutoff_derivative, stratmann_cutoff_derivative
   public :: becke_cutoff_second_derivative, stratmann_cutoff_second_derivative
   public :: partition_scheme_name

   ! Cutoff profile
   integer, parameter :: PARTITION_BECKE = 1     !! Three iterations of p(x), smooth everywhere
   integer, parameter :: PARTITION_STRATMANN = 2  !! Degree-5 with a hard cutoff, screenable

   ! Atomic size adjustment
   integer, parameter :: ADJUST_NONE = 0      !! All atoms the same size
   integer, parameter :: ADJUST_BECKE = 1     !! Bragg radii
   integer, parameter :: ADJUST_TREUTLER = 2  !! Square roots of the Bragg radii

   real(dp), parameter :: STRATMANN_A = 0.64_dp
      !! Stratmann's cutoff parameter, eq. 14

   real(dp), parameter :: MAX_ADJUST = 0.5_dp
      !! Becke's clamp on the size-adjustment shift

   real(dp), parameter, public :: TINY_RADIUS = 1.0e-200_dp
      !! Guards the ratio when a radius is zero (a ghost atom)

   integer, parameter :: N_DIM = 3

   integer, parameter :: SCREEN_NONE = 0       !! Every pair, the reference
   integer, parameter :: SCREEN_STRATMANN = 1  !! Exact: skip what is exactly zero
   integer, parameter :: SCREEN_BOUNDED = 2    !! Becke: drop what is provably negligible

   integer, parameter :: MIN_STRATMANN_ATOMS = 12
      !! Below this a point costs less to evaluate in full than to classify
   integer, parameter :: MIN_BOUNDED_ATOMS = 48
      !! The Becke bound needs a survivor's whole product on top of the bounds, and
      !! measured slower than the plain loop at 30 atoms (1.7 s against 0.8 s
      !! serial on ten waters); it is faster at 96
   integer, parameter :: N_NEAR = 6
      !! Atoms nearest a point used to bound the Becke cell functions
   integer, parameter :: POINT_CHUNK = 64
      !! Grid points handed to a thread at a time

   real(dp), parameter :: SCREEN_TOL = 1.0e-18_dp
      !! Becke screening drops an atom whose cell function is provably below
      !! SCREEN_TOL times the owner's. Each weight then moves by a relative
      !! (n_atoms - 1) * SCREEN_TOL at most: 2e-15 at 2000 atoms, inside the 1e-14
      !! this scheme promises up to that size with a factor of five to spare.

contains

   pure function partition_scheme_name(scheme) result(name)
      !! Human-readable scheme name, for logs and error messages
      integer, intent(in) :: scheme
      character(len=:), allocatable :: name

      select case (scheme)
      case (PARTITION_BECKE)
         name = "Becke"
      case (PARTITION_STRATMANN)
         name = "Stratmann"
      case default
         name = "unknown"
      end select
   end function partition_scheme_name

   subroutine becke_partition_weights(points, atom_coords, atomic_numbers, owner, &
                                      scheme, adjust, weights, error, screening)
      !! Partition weight at each grid point, for the atom that owns it
      !!
      !! The returned weight is w_owner(r): the fraction of the integrand at
      !! that point assigned to the atom whose atomic grid produced it. Multiply
      !! it by the point's radial and angular weight to get its contribution.
      !!
      !! Threaded over points. Every point is independent and writes only its own
      !! `weights(k)`, so the result does not depend on the thread count.
      !!
      !! Screening (on by default) leaves the weights unchanged under
      !! `PARTITION_STRATMANN` (bitwise) and moves each one by a relative
      !! `n_atoms * SCREEN_TOL` at most under `PARTITION_BECKE`; see
      !! `cell_stratmann` and `cell_bounded`.
      real(dp), intent(in) :: points(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: atom_coords(:, :)  !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)   !! Z per atom, for the radii
      integer, intent(in) :: owner(:)            !! Atom each point belongs to
      integer, intent(in) :: scheme              !! PARTITION_BECKE or PARTITION_STRATMANN
      integer, intent(in) :: adjust              !! ADJUST_NONE, _BECKE or _TREUTLER
      real(dp), intent(out) :: weights(:)        !! (n_points)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: screening
         !! `.false.` evaluates every atom pair at every point, the reference the
         !! screened paths are checked against. Default `.true.`

      real(dp), allocatable :: shift(:, :), inv_distance(:, :), cell(:), atom_r(:)
      integer, allocatable :: near(:)
      integer :: n_atoms, n_points, i, j, k, mode
      real(dp) :: total

      n_atoms = size(atom_coords, 2)
      n_points = size(points, 2)

      if (size(atomic_numbers) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "partition: atomic_numbers does not match atom_coords")
         return
      end if
      if (size(owner) /= n_points .or. size(weights) /= n_points) then
         call error%set(ERROR_VALIDATION, "partition: owner and weights must match the points")
         return
      end if
      if (scheme /= PARTITION_BECKE .and. scheme /= PARTITION_STRATMANN) then
         call error%set(ERROR_VALIDATION, "partition: unknown scheme")
         return
      end if
      if (any(owner < 1) .or. any(owner > n_atoms)) then
         call error%set(ERROR_VALIDATION, "partition: owner index outside the atom list")
         return
      end if

      ! A single atom owns everything, and the pair loop below would leave the
      ! product empty. Short-circuit rather than special-case inside the loop.
      if (n_atoms == 1) then
         weights = 1.0_dp
         return
      end if

      mode = SCREEN_NONE
      if (scheme == PARTITION_STRATMANN .and. n_atoms >= MIN_STRATMANN_ATOMS) then
         mode = SCREEN_STRATMANN
      else if (scheme == PARTITION_BECKE .and. n_atoms >= MIN_BOUNDED_ATOMS) then
         mode = SCREEN_BOUNDED
      end if
      if (present(screening)) then
         if (.not. screening) mode = SCREEN_NONE
      end if

      allocate (shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms))

      call size_adjustment(atomic_numbers, adjust, shift)
      ! Only the upper triangle is ever the shift of a pair (lo, hi); mirroring
      ! it makes the matrix symmetric for `pair_factor`'s column access.
      do i = 1, n_atoms
         do j = 1, i - 1
            shift(i, j) = shift(j, i)
         end do
      end do

      ! 1/R_AB, precomputed: it is the same for every point.
      inv_distance = 0.0_dp
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) then
               inv_distance(i, j) = 1.0_dp/norm2(atom_coords(:, i) - atom_coords(:, j))
            end if
         end do
      end do

      ! Scratch is allocated inside the region so each thread owns its own;
      ! `shift` and `inv_distance` depend only on the nuclei and stay shared.
      ! Dynamic, because a screened point costs between one atom's product and
      ! a handful, depending on how far it is from the nearest nucleus.
      !$omp parallel default(none) &
      !$omp    shared(points, atom_coords, owner, scheme, weights, shift, &
      !$omp           inv_distance, n_points, n_atoms, mode) &
      !$omp    private(k, i, atom_r, cell, near, total)
      allocate (atom_r(n_atoms), cell(n_atoms), near(min(N_NEAR, n_atoms)))

      !$omp do schedule(dynamic, POINT_CHUNK)
      do k = 1, n_points
         do i = 1, n_atoms
            atom_r(i) = norm2(points(:, k) - atom_coords(:, i))
         end do

         select case (mode)
         case (SCREEN_STRATMANN)
            call cell_stratmann(atom_r, inv_distance, shift, cell)
         case (SCREEN_BOUNDED)
            call cell_bounded(owner(k), atom_r, inv_distance, shift, near, cell)
         case default
            call cell_full(scheme, atom_r, inv_distance, shift, cell)
         end select

         total = sum(cell)
         if (total > 0.0_dp) then
            weights(k) = cell(owner(k))/total
         else
            ! Every cell function underflowed, which happens only absurdly far
            ! from the molecule where the integrand is zero anyway.
            weights(k) = 0.0_dp
         end if
      end do
      !$omp end do
      !$omp end parallel
   end subroutine becke_partition_weights

   pure subroutine cell_full(scheme, atom_r, inv_distance, shift, cell)
      !! Unnormalised cell function of every atom at one point, every pair visited
      integer, intent(in) :: scheme
      real(dp), intent(in) :: atom_r(:), inv_distance(:, :), shift(:, :)
      real(dp), intent(out) :: cell(:)

      integer :: i, j
      real(dp) :: s

      cell = 1.0_dp
      do i = 1, size(cell)
         do j = i + 1, size(cell)
            s = cutoff(scheme, pair_nu(i, j, atom_r, inv_distance, shift))
            ! s is the fraction going to i; the complement goes to j. Doing
            ! both here keeps the pair symmetric by construction, so the
            ! weights cannot fail to sum to 1 through an asymmetry in nu.
            cell(i) = cell(i)*s
            cell(j) = cell(j)*(1.0_dp - s)
         end do
      end do
   end subroutine cell_full

   pure subroutine cell_stratmann(atom_r, inv_distance, shift, cell)
      !! Cell functions at one point, skipping what is exactly zero
      !!
      !! Exact: the result is bitwise that of `cell_full`. Let m be the atom
      !! nearest the point. If the cutoff of the pair (j, m) is flat at 0 for j,
      !! then cell(j) carries a factor of exactly zero and is exactly zero, so j
      !! is skipped. Every atom that survives is multiplied out over *all* its
      !! partners, in the same order `cell_full` uses (partner index ascending),
      !! because a dead partner can still contribute a fractional factor to a
      !! live atom. Cost is O(n_atoms) to classify plus O(n_atoms) per survivor.
      real(dp), intent(in) :: atom_r(:), inv_distance(:, :), shift(:, :)
      real(dp), intent(out) :: cell(:)

      integer :: n, m, j, k
      real(dp) :: nu, p

      n = size(cell)
      m = minloc(atom_r, dim=1)
      cell = 0.0_dp
      do j = 1, n
         if (j < m) then
            nu = pair_nu(j, m, atom_r, inv_distance, shift)
            if (nu >= STRATMANN_A) cycle
         else if (j > m) then
            nu = pair_nu(m, j, atom_r, inv_distance, shift)
            if (nu <= -STRATMANN_A) cycle
         end if
         p = 1.0_dp
         do k = 1, n
            if (k /= j) p = p*pair_factor(PARTITION_STRATMANN, j, k, atom_r, inv_distance, shift)
         end do
         cell(j) = p
      end do
   end subroutine cell_stratmann

   pure subroutine cell_bounded(owner, atom_r, inv_distance, shift, near, cell)
      !! Cell functions at one point, dropping atoms whose cell is negligible
      !!
      !! Every factor of P_i = prod_{j/=i} s(nu_ij) lies in [0, 1], so a product
      !! over any subset of the partners bounds it from above. Taking the
      !! `N_NEAR` atoms nearest the point as that subset (they are the ones that
      !! make the factors small), an atom i is dropped when the bound
      !! b_i <= SCREEN_TOL * P_owner, and the owner's own P is always computed
      !! in full. Each dropped atom then removes at most SCREEN_TOL * P_owner
      !! from sum_i P_i, which is itself at least P_owner, so the weight
      !! w = P_owner / sum_i P_i changes by a relative amount of at most
      !! (n_atoms - 1) * SCREEN_TOL. The same test is applied to the running
      !! product while a survivor is multiplied out, which is again a subset.
      !! Survivors are multiplied in the order `cell_full` uses, so they are
      !! bitwise its values; only the dropped mass differs.
      integer, intent(in) :: owner
      real(dp), intent(in) :: atom_r(:), inv_distance(:, :), shift(:, :)
      integer, intent(inout) :: near(:)
      real(dp), intent(out) :: cell(:)

      integer :: n, kn, i, j, q, slot
      real(dp) :: p_owner, thresh, b
      logical :: dropped

      n = size(cell)
      kn = size(near)

      ! The kn smallest distances, kept sorted, by insertion.
      near = 0
      do i = 1, n
         if (near(kn) /= 0) then
            if (atom_r(i) >= atom_r(near(kn))) cycle
         end if
         slot = kn
         do while (slot > 1)
            if (near(slot - 1) /= 0) then
               if (atom_r(near(slot - 1)) <= atom_r(i)) exit
            end if
            near(slot) = near(slot - 1)
            slot = slot - 1
         end do
         near(slot) = i
      end do

      p_owner = 1.0_dp
      do j = 1, n
         if (j /= owner) then
            p_owner = p_owner*pair_factor(PARTITION_BECKE, owner, j, atom_r, inv_distance, shift)
         end if
      end do
      thresh = SCREEN_TOL*p_owner

      cell = 0.0_dp
      cell(owner) = p_owner
      do i = 1, n
         if (i == owner) cycle
         dropped = .false.
         b = 1.0_dp
         do q = 1, kn
            j = near(q)
            if (j == i) cycle
            b = b*pair_factor(PARTITION_BECKE, i, j, atom_r, inv_distance, shift)
            if (b <= thresh) then
               dropped = .true.
               exit
            end if
         end do
         if (dropped) cycle

         b = 1.0_dp
         do j = 1, n
            if (j == i) cycle
            b = b*pair_factor(PARTITION_BECKE, i, j, atom_r, inv_distance, shift)
            if (b <= thresh) then
               dropped = .true.
               exit
            end if
         end do
         if (.not. dropped) cell(i) = b
      end do
   end subroutine cell_bounded

   pure function pair_nu(lo, hi, atom_r, inv_distance, shift) result(nu)
      !! The size-adjusted elliptical coordinate of the pair (lo, hi), lo < hi
      integer, intent(in) :: lo, hi
      real(dp), intent(in) :: atom_r(:), inv_distance(:, :), shift(:, :)
      real(dp) :: nu
      real(dp) :: mu

      mu = (atom_r(lo) - atom_r(hi))*inv_distance(lo, hi)
      nu = mu + shift(lo, hi)*(1.0_dp - mu*mu)
   end function pair_nu

   pure function cutoff(scheme, nu) result(s)
      !! The cutoff of `scheme` at `nu`
      integer, intent(in) :: scheme
      real(dp), intent(in) :: nu
      real(dp) :: s

      if (scheme == PARTITION_BECKE) then
         s = becke_cutoff(nu)
      else
         s = stratmann_cutoff(nu)
      end if
   end function cutoff

   pure function pair_factor(scheme, i, j, atom_r, inv_distance, shift) result(f)
      !! The factor atom i's cell function takes from partner j
      !!
      !! s(nu_ij) when i < j and 1 - s(nu_ji) when i > j, the same two numbers
      !! `cell_full` multiplies in, so a product of these is bitwise its cell.
      !! `inv_distance` and `shift` are symmetric here (the triangle `pair_nu`
      !! reads is mirrored by the caller), which lets a loop over j at fixed i
      !! walk down a column instead of across a row.
      integer, intent(in) :: scheme, i, j
      real(dp), intent(in) :: atom_r(:), inv_distance(:, :), shift(:, :)
      real(dp) :: f
      real(dp) :: mu, nu

      if (i < j) then
         mu = (atom_r(i) - atom_r(j))*inv_distance(j, i)
         nu = mu + shift(j, i)*(1.0_dp - mu*mu)
         f = cutoff(scheme, nu)
      else
         mu = (atom_r(j) - atom_r(i))*inv_distance(j, i)
         nu = mu + shift(j, i)*(1.0_dp - mu*mu)
         f = 1.0_dp - cutoff(scheme, nu)
      end if
   end function pair_factor

   subroutine becke_partition_derivatives(points, atom_coords, atomic_numbers, owner, &
                                          scheme, adjust, dweights, error, screening)
      !! dP_k/dR_A at fixed grid points, as (3, n_atoms, n_points)
      !!
      !! The partition weight depends on every nuclear position, so moving any
      !! atom changes the weight of every point, including points owned by other
      !! atoms. That is the "grid response" term of a DFT gradient.
      !!
      !! **The grid points are held fixed here.** Points also move with the atom
      !! that owns them; that is a separate contribution the caller adds.
      !!
      !! The derivative is carried through the cell product alongside the value
      !! rather than as a logarithmic derivative summed at the end: the log form
      !! divides by each cutoff factor, and those reach zero far from the
      !! molecule where the derivative does not, so the quotient is 0/0 exactly
      !! where a guard would have to guess.
      !!
      !! Screening (on by default) is exact and applies to `PARTITION_STRATMANN`
      !! only. An atom whose cutoff against the nearest atom is flat at zero has a
      !! cell function *and* a cell derivative of exactly zero, since every term of
      !! the product rule contains either that zero factor or that flat
      !! derivative, so it is skipped, as is any partner whose factor is exactly 1
      !! with an exactly zero derivative. The result is bitwise the unscreened one.
      !! `PARTITION_BECKE` has no exact zeros to skip, and a bound for it would
      !! need a bound on the derivative as well as on the value; it is evaluated
      !! in full.
      real(dp), intent(in) :: points(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: atom_coords(:, :)  !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)   !! Z per atom, for the radii
      integer, intent(in) :: owner(:)            !! Atom each point belongs to
      integer, intent(in) :: scheme
      integer, intent(in) :: adjust
      real(dp), intent(out) :: dweights(:, :, :)  !! (3, n_atoms, n_points)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: screening
         !! `.false.` keeps every atom as a candidate. Default `.true.`

      real(dp), allocatable :: shift(:, :), inv_distance(:, :)
      real(dp), allocatable :: atom_r(:), unit_vec(:, :)
      real(dp), allocatable :: cell(:), dcell(:, :, :), grow(:, :, :)
      integer, allocatable :: cand(:), slot(:)
      real(dp) :: dtotal(3), rij(3), ds_lo(3), ds_hi(3), dmu_lo(3), dmu_hi(3)
      real(dp) :: df_i(3), df_k(3)
      real(dp) :: mu, nu, s, f, total, dnu_dmu, ds_dnu, invd, ci
      integer :: n_atoms, n_points, i, j, k, a, c, nc, m, lo, hi
      logical :: screen, dead

      n_atoms = size(atom_coords, 2)
      n_points = size(points, 2)
      dweights = 0.0_dp

      if (size(atomic_numbers) /= n_atoms) then
         call error%set(ERROR_VALIDATION, &
                        "partition derivatives: atomic_numbers does not match atom_coords")
         return
      end if
      if (size(owner) /= n_points .or. size(dweights, 3) /= n_points &
          .or. size(dweights, 2) /= n_atoms) then
         call error%set(ERROR_VALIDATION, &
                        "partition derivatives: owner and dweights must match the points and atoms")
         return
      end if
      if (scheme /= PARTITION_BECKE .and. scheme /= PARTITION_STRATMANN) then
         call error%set(ERROR_VALIDATION, "partition derivatives: unknown scheme")
         return
      end if
      if (any(owner < 1) .or. any(owner > n_atoms)) then
         call error%set(ERROR_VALIDATION, "partition derivatives: owner index outside the atom list")
         return
      end if

      ! One atom owns everything and its weight is identically one, so it does
      ! not change when the atom moves.
      if (n_atoms == 1) return

      screen = scheme == PARTITION_STRATMANN
      if (present(screening)) screen = screen .and. screening

      allocate (shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms))

      call size_adjustment(atomic_numbers, adjust, shift)

      inv_distance = 0.0_dp
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) then
               inv_distance(i, j) = 1.0_dp/norm2(atom_coords(:, i) - atom_coords(:, j))
            end if
         end do
      end do

      ! One thread per grid point. Every iteration writes its own
      ! `dweights(:, :, k)`, so nothing is reduced; the per-point scratch is
      ! allocated inside the region so each thread gets its own, while `shift`
      ! and `inv_distance` depend only on the nuclei and stay shared. Dynamic
      ! scheduling because a screened point's cost depends on where it sits.
      !
      ! `dcell(:, :, c)` is the derivative of candidate c's cell function with
      ! respect to every nucleus. It is sized to the candidate count and grown on
      ! demand, so a screened run does not carry an (n_atoms, n_atoms) block per
      ! thread.
      !$omp parallel default(none) &
      !$omp    shared(points, atom_coords, owner, scheme, dweights, shift, &
      !$omp           inv_distance, n_points, n_atoms, screen) &
      !$omp    private(k, i, j, a, c, nc, m, lo, hi, dead, atom_r, unit_vec, cell, dcell, grow, &
      !$omp            cand, slot, dtotal, rij, ds_lo, ds_hi, dmu_lo, dmu_hi, df_i, df_k, &
      !$omp            mu, nu, s, f, total, dnu_dmu, ds_dnu, invd, ci)
      allocate (atom_r(n_atoms), unit_vec(3, n_atoms), cell(n_atoms))
      allocate (cand(n_atoms), slot(n_atoms))
      allocate (dcell(3, n_atoms, min(n_atoms, 8)))

      !$omp do schedule(dynamic, POINT_CHUNK)
      do k = 1, n_points
         do i = 1, n_atoms
            atom_r(i) = norm2(points(:, k) - atom_coords(:, i))
            if (atom_r(i) > 0.0_dp) then
               unit_vec(:, i) = (points(:, k) - atom_coords(:, i))/atom_r(i)
            else
               ! The point sits on the nucleus. |r - R| is not differentiable
               ! there, and the quadrature never places a point exactly on one,
               ! so this is a guard rather than a case.
               unit_vec(:, i) = 0.0_dp
            end if
         end do

         ! The atoms whose cell function can be non-zero, in ascending order.
         nc = 0
         slot = 0
         m = minloc(atom_r, dim=1)
         do i = 1, n_atoms
            dead = .false.
            if (screen) then
               if (i < m) then
                  dead = pair_nu(i, m, atom_r, inv_distance, shift) >= STRATMANN_A
               else if (i > m) then
                  dead = pair_nu(m, i, atom_r, inv_distance, shift) <= -STRATMANN_A
               end if
            end if
            if (.not. dead) then
               nc = nc + 1
               cand(nc) = i
               slot(i) = nc
            end if
         end do

         if (nc > size(dcell, 3)) then
            allocate (grow(3, n_atoms, min(n_atoms, max(nc, 2*size(dcell, 3)))))
            call move_alloc(grow, dcell)
         end if

         cell = 0.0_dp
         do c = 1, nc
            i = cand(c)
            ci = 1.0_dp
            dcell(:, :, c) = 0.0_dp

            ! Partners in ascending order, which is the order the pairs reach
            ! atom i in a double loop over (lo, hi): the value and derivative
            ! then accumulate in the same sequence as an unscreened run.
            do j = 1, n_atoms
               if (j == i) cycle
               lo = min(i, j)
               hi = max(i, j)
               invd = inv_distance(lo, hi)

               mu = (atom_r(lo) - atom_r(hi))*invd
               nu = mu + shift(lo, hi)*(1.0_dp - mu*mu)

               if (scheme == PARTITION_BECKE) then
                  s = becke_cutoff(nu)
                  ds_dnu = becke_cutoff_derivative(nu)
               else
                  s = stratmann_cutoff(nu)
                  ds_dnu = stratmann_cutoff_derivative(nu)
               end if

               if (i < j) then
                  f = s
               else
                  f = 1.0_dp - s
               end if
               if (f == 1.0_dp .and. ds_dnu == 0.0_dp) cycle

               rij = atom_coords(:, lo) - atom_coords(:, hi)
               dnu_dmu = 1.0_dp - 2.0_dp*shift(lo, hi)*mu

               ! mu depends on R_lo and R_hi only: through the distances from the
               ! point, and through the internuclear distance in the denominator.
               dmu_lo = -unit_vec(:, lo)*invd - (atom_r(lo) - atom_r(hi))*invd**3*rij
               dmu_hi = unit_vec(:, hi)*invd + (atom_r(lo) - atom_r(hi))*invd**3*rij

               ds_lo = ds_dnu*dnu_dmu*dmu_lo
               ds_hi = ds_dnu*dnu_dmu*dmu_hi

               if (i < j) then
                  df_i = ds_lo
                  df_k = ds_hi
               else
                  df_i = -ds_hi
                  df_k = -ds_lo
               end if

               ! Product rule, value and derivative together. Order matters:
               ! the derivative uses the cell value from before this factor.
               dcell(:, :, c) = dcell(:, :, c)*f
               dcell(:, i, c) = dcell(:, i, c) + ci*df_i
               dcell(:, j, c) = dcell(:, j, c) + ci*df_k
               ci = ci*f
            end do
            cell(i) = ci
         end do

         total = sum(cell)
         if (total <= 0.0_dp) cycle
         ! A skipped owner has a cell function and derivative of exactly zero.
         if (slot(owner(k)) == 0) cycle

         do a = 1, n_atoms
            dtotal = 0.0_dp
            do c = 1, nc
               dtotal = dtotal + dcell(:, a, c)
            end do
            dweights(:, a, k) = (dcell(:, a, slot(owner(k)))*total &
                                 - cell(owner(k))*dtotal)/(total*total)
         end do
      end do
      !$omp end do
      !$omp end parallel
   end subroutine becke_partition_derivatives

   pure subroutine size_adjustment(atomic_numbers, adjust, shift)
      !! Pairwise shift applied to mu so unequal atoms get unequal cells
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: adjust
      real(dp), intent(out) :: shift(:, :)

      real(dp), allocatable :: radius(:)
      integer :: n_atoms, i, j
      real(dp) :: chi, a

      n_atoms = size(atomic_numbers)
      allocate (radius(n_atoms))

      select case (adjust)
      case (ADJUST_TREUTLER)
         do i = 1, n_atoms
            radius(i) = sqrt(bragg_radius(atomic_numbers(i))) + TINY_RADIUS
         end do
      case (ADJUST_BECKE)
         do i = 1, n_atoms
            radius(i) = bragg_radius(atomic_numbers(i)) + TINY_RADIUS
         end do
      case default
         radius = 1.0_dp
      end select

      shift = 0.0_dp
      if (adjust == ADJUST_NONE) return

      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i == j) cycle
            chi = radius(i)/radius(j)
            ! a = u/(u^2-1) with u = (chi-1)/(chi+1), which reduces to this.
            a = 0.25_dp*(1.0_dp/chi - chi)
            shift(i, j) = max(-MAX_ADJUST, min(MAX_ADJUST, a))
         end do
      end do
   end subroutine size_adjustment

   pure function becke_cutoff(nu) result(s)
      !! Becke's smooth cutoff: three iterations of p(x) = x(3 - x^2)/2
      real(dp), intent(in) :: nu
      real(dp) :: s
      real(dp) :: f

      f = nu
      f = 0.5_dp*f*(3.0_dp - f*f)
      f = 0.5_dp*f*(3.0_dp - f*f)
      f = 0.5_dp*f*(3.0_dp - f*f)
      s = 0.5_dp*(1.0_dp - f)
   end function becke_cutoff

   pure function stratmann_cutoff(mu) result(s)
      !! Stratmann's cutoff, exactly 1 or 0 outside |mu| >= a
      real(dp), intent(in) :: mu
      real(dp) :: s
      real(dp) :: z, z2

      if (mu <= -STRATMANN_A) then
         s = 1.0_dp
      else if (mu >= STRATMANN_A) then
         s = 0.0_dp
      else
         z = mu/STRATMANN_A
         z2 = z*z
         z = (1.0_dp/16.0_dp)*(z*(35.0_dp + z2*(-35.0_dp + z2*(21.0_dp - 5.0_dp*z2))))
         s = 0.5_dp*(1.0_dp - z)
      end if
   end function stratmann_cutoff

   pure function becke_cutoff_derivative(nu) result(ds)
      !! d/dnu of `becke_cutoff`
      !!
      !! Each iteration is p(x) = x(3 - x^2)/2 with p'(x) = 3(1 - x^2)/2, so the
      !! chain rule over the three of them is the product of p' at the value
      !! going into each.
      real(dp), intent(in) :: nu
      real(dp) :: ds
      real(dp) :: f0, f1, f2

      f0 = nu
      f1 = 0.5_dp*f0*(3.0_dp - f0*f0)
      f2 = 0.5_dp*f1*(3.0_dp - f1*f1)

      ds = -0.5_dp*(1.5_dp*(1.0_dp - f2*f2)) &
           *(1.5_dp*(1.0_dp - f1*f1)) &
           *(1.5_dp*(1.0_dp - f0*f0))
   end function becke_cutoff_derivative

   pure function stratmann_cutoff_derivative(mu) result(ds)
      !! d/dmu of `stratmann_cutoff`, zero where the cutoff is flat
      real(dp), intent(in) :: mu
      real(dp) :: ds
      real(dp) :: z, z2

      if (mu <= -STRATMANN_A .or. mu >= STRATMANN_A) then
         ds = 0.0_dp
      else
         z = mu/STRATMANN_A
         z2 = z*z
         ds = -0.5_dp*(1.0_dp/16.0_dp) &
              *(35.0_dp + z2*(-105.0_dp + z2*(105.0_dp - 35.0_dp*z2)))/STRATMANN_A
      end if
   end function stratmann_cutoff_derivative

   pure function becke_cutoff_second_derivative(nu) result(d2s)
      !! d2/dnu2 of `becke_cutoff`
      !!
      !! Three nested copies of p(x) = x(3 - x^2)/2, so with
      !! f0 = nu, f1 = p(f0), f2 = p(f1) and s = (1 - p(f2))/2 the first
      !! derivative is the product of the three p' factors and the second
      !! follows from differentiating that product:
      !!
      !!   s'  = -1/2 . p'(f2) p'(f1) p'(f0)
      !!   s'' = -1/2 [ p''(f2) p'(f1)^2 p'(f0)^2
      !!                + p'(f2) p''(f1) p'(f0)^2
      !!                + p'(f2) p'(f1) p''(f0) ]
      !!
      !! using p'(x) = 3(1 - x^2)/2 and p''(x) = -3x. Each term carries the
      !! chain factors for the layers *below* the one being differentiated
      !! twice, which is where the squares come from.
      real(dp), intent(in) :: nu
      real(dp) :: d2s
      real(dp) :: f0, f1, f2, p0, p1, p2, q0, q1, q2

      f0 = nu
      f1 = 0.5_dp*f0*(3.0_dp - f0*f0)
      f2 = 0.5_dp*f1*(3.0_dp - f1*f1)

      p0 = 1.5_dp*(1.0_dp - f0*f0)
      p1 = 1.5_dp*(1.0_dp - f1*f1)
      p2 = 1.5_dp*(1.0_dp - f2*f2)

      q0 = -3.0_dp*f0
      q1 = -3.0_dp*f1
      q2 = -3.0_dp*f2

      d2s = -0.5_dp*(q2*p1*p1*p0*p0 + p2*q1*p0*p0 + p2*p1*q0)
   end function becke_cutoff_second_derivative

   pure function stratmann_cutoff_second_derivative(mu) result(d2s)
      !! d2/dmu2 of `stratmann_cutoff`, zero where the cutoff is flat
      !!
      !! The polynomial is g(z) = (35z - 35z^3 + 21z^5 - 5z^7)/16 with
      !! s = (1 - g)/2, so s'' = -g''/2 scaled by 1/a^2 for the two chain
      !! factors of z = mu/a.
      real(dp), intent(in) :: mu
      real(dp) :: d2s
      real(dp) :: z, z2

      if (mu <= -STRATMANN_A .or. mu >= STRATMANN_A) then
         d2s = 0.0_dp
      else
         z = mu/STRATMANN_A
         z2 = z*z
         d2s = -0.5_dp*(1.0_dp/16.0_dp) &
               *z*(-210.0_dp + z2*(420.0_dp - 210.0_dp*z2)) &
               /(STRATMANN_A*STRATMANN_A)
      end if
   end function stratmann_cutoff_second_derivative

end module mqc_dft_partition
