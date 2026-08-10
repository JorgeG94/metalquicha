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
   !! which is smooth everywhere and non-zero everywhere -- so every atom
   !! contributes at every point, and the cost is unavoidably natoms^2 per
   !! point. Stratmann's [CPL 257, 213 (1996)] is a degree-5 polynomial that
   !! reaches exactly +-1 at |mu| >= a, so distant pairs contribute exactly 0 or
   !! 1 and can be skipped. That is what makes large systems affordable.
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
   public :: partition_scheme_name

   !> Cutoff profile
   integer, parameter :: PARTITION_BECKE = 1     !! Three iterations of p(x), smooth everywhere
   integer, parameter :: PARTITION_STRATMANN = 2  !! Degree-5 with a hard cutoff, screenable

   !> Atomic size adjustment
   integer, parameter :: ADJUST_NONE = 0      !! All atoms the same size
   integer, parameter :: ADJUST_BECKE = 1     !! Bragg radii
   integer, parameter :: ADJUST_TREUTLER = 2  !! Square roots of the Bragg radii

   !> Stratmann's cutoff parameter, eq. 14
   real(dp), parameter :: STRATMANN_A = 0.64_dp

   !> Becke's clamp on the size-adjustment shift
   real(dp), parameter :: MAX_ADJUST = 0.5_dp

   !> Guards the ratio when a radius is zero (a ghost atom)
   real(dp), parameter :: TINY_RADIUS = 1.0e-200_dp

   integer, parameter :: N_DIM = 3

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
                                      scheme, adjust, weights, error)
      !! Partition weight at each grid point, for the atom that owns it
      !!
      !! The returned weight is w_owner(r): the fraction of the integrand at
      !! that point assigned to the atom whose atomic grid produced it. Multiply
      !! it by the point's radial and angular weight to get its contribution.
      real(dp), intent(in) :: points(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: atom_coords(:, :)  !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)   !! Z per atom, for the radii
      integer, intent(in) :: owner(:)            !! Atom each point belongs to
      integer, intent(in) :: scheme              !! PARTITION_BECKE or PARTITION_STRATMANN
      integer, intent(in) :: adjust              !! ADJUST_NONE, _BECKE or _TREUTLER
      real(dp), intent(out) :: weights(:)        !! (n_points)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: shift(:, :), inv_distance(:, :), cell(:), atom_r(:)
      integer :: n_atoms, n_points, i, j, k
      real(dp) :: mu, nu, s, total

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

      allocate (shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms))
      allocate (cell(n_atoms), atom_r(n_atoms))

      call size_adjustment(atomic_numbers, adjust, shift)

      ! 1/R_AB, precomputed: it is the same for every point.
      inv_distance = 0.0_dp
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) then
               inv_distance(i, j) = 1.0_dp/norm2(atom_coords(:, i) - atom_coords(:, j))
            end if
         end do
      end do

      do k = 1, n_points
         do i = 1, n_atoms
            atom_r(i) = norm2(points(:, k) - atom_coords(:, i))
         end do

         cell = 1.0_dp
         do i = 1, n_atoms
            do j = i + 1, n_atoms
               mu = (atom_r(i) - atom_r(j))*inv_distance(i, j)
               nu = mu + shift(i, j)*(1.0_dp - mu*mu)

               if (scheme == PARTITION_BECKE) then
                  s = becke_cutoff(nu)
               else
                  s = stratmann_cutoff(nu)
               end if

               ! s is the fraction going to i; the complement goes to j. Doing
               ! both here halves the work and keeps the pair symmetric by
               ! construction, so the weights cannot fail to sum to 1 through
               ! an asymmetry in nu.
               cell(i) = cell(i)*s
               cell(j) = cell(j)*(1.0_dp - s)
            end do
         end do

         total = sum(cell)
         if (total > 0.0_dp) then
            weights(k) = cell(owner(k))/total
         else
            ! Every cell function underflowed, which happens only absurdly far
            ! from the molecule where the integrand is zero anyway.
            weights(k) = 0.0_dp
         end if
      end do
   end subroutine becke_partition_weights

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

end module mqc_dft_partition
