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
   public :: becke_partition_derivatives
   public :: becke_cutoff_derivative, stratmann_cutoff_derivative
   public :: becke_cutoff_second_derivative, stratmann_cutoff_second_derivative
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
   real(dp), parameter, public :: TINY_RADIUS = 1.0e-200_dp

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

   subroutine becke_partition_derivatives(points, atom_coords, atomic_numbers, owner, &
                                          scheme, adjust, dweights, error)
      !! dP_k/dR_A at fixed grid points, as (3, n_atoms, n_points)
      !!
      !! The partition weight depends on every nuclear position, so moving any
      !! atom changes the weight of every point -- including points owned by
      !! other atoms. That is the whole content of the "grid response" term in
      !! a DFT gradient, and omitting it leaves a gradient that looks entirely
      !! reasonable and disagrees with finite differences at about 1e-4.
      !!
      !! **The grid points are held fixed here.** Points move with the atom that
      !! owns them, but that motion is a separate contribution which the caller
      !! adds; mixing the two into one routine would make neither checkable.
      !!
      !! The cell functions are differentiated by carrying the derivative
      !! through the product alongside the value, rather than as a logarithmic
      !! derivative summed at the end. The log form needs a division by each
      !! cutoff factor, and those legitimately reach zero far from the molecule
      !! -- where the value underflows to zero but the derivative does not, so
      !! the quotient is 0/0 exactly where a guard would have to guess. Carrying
      !! both costs an extra factor of the atom count and cannot divide by
      !! anything.
      real(dp), intent(in) :: points(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: atom_coords(:, :)  !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)   !! Z per atom, for the radii
      integer, intent(in) :: owner(:)            !! Atom each point belongs to
      integer, intent(in) :: scheme
      integer, intent(in) :: adjust
      real(dp), intent(out) :: dweights(:, :, :)  !! (3, n_atoms, n_points)
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: shift(:, :), inv_distance(:, :)
      real(dp), allocatable :: atom_r(:), unit_vec(:, :)
      real(dp), allocatable :: cell(:), dcell(:, :, :)
      real(dp) :: dtotal(3), rij(3), ds_i(3), ds_j(3), dmu_i(3), dmu_j(3)
      real(dp) :: mu, nu, s, total, dnu_dmu, ds_dnu, invd
      integer :: n_atoms, n_points, i, j, k, a

      n_atoms = size(atom_coords, 2)
      n_points = size(points, 2)
      dweights = 0.0_dp

      if (size(atomic_numbers) /= n_atoms) then
         call error%set(ERROR_VALIDATION, &
                        "partition derivatives: atomic_numbers does not match atom_coords")
         return
      end if
      if (scheme /= PARTITION_BECKE .and. scheme /= PARTITION_STRATMANN) then
         call error%set(ERROR_VALIDATION, "partition derivatives: unknown scheme")
         return
      end if

      ! One atom owns everything and its weight is identically one, so it does
      ! not change when the atom moves.
      if (n_atoms == 1) return

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

      ! One thread per grid point.
      !
      ! Every iteration reads the geometry and writes its own `dweights(:, :, k)`,
      ! so the points are independent and nothing is reduced. What is not
      ! independent is the scratch -- `cell` and `dcell` are rebuilt per point --
      ! so those move inside the region for each thread to allocate its own.
      ! `shift` and `inv_distance` depend only on the nuclei and stay shared.
      !
      ! This is the grid-response term of a density functional gradient, and it
      ! was 100 s of CPU on twenty waters while this file had no OpenMP in it at
      ! all. `dcell` is the only sizeable copy, three by the atom count squared:
      ! 86 kB for twenty waters, so sixteen threads cost under two megabytes.
      !
      ! `schedule(static)`, unlike the quadrature loops: every point costs the
      ! same `n_atoms` squared wherever it sits, so there is nothing for dynamic
      ! scheduling to balance and its bookkeeping would be pure overhead.
      !$omp parallel default(none) &
      !$omp    shared(points, atom_coords, owner, scheme, dweights, shift, &
      !$omp           inv_distance, n_points, n_atoms) &
      !$omp    private(k, i, j, a, atom_r, unit_vec, cell, dcell, dtotal, rij, &
      !$omp            ds_i, ds_j, dmu_i, dmu_j, mu, nu, s, total, dnu_dmu, &
      !$omp            ds_dnu, invd)
      allocate (atom_r(n_atoms), unit_vec(3, n_atoms))
      allocate (cell(n_atoms), dcell(3, n_atoms, n_atoms))

      !$omp do schedule(static)
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

         cell = 1.0_dp
         dcell = 0.0_dp

         do i = 1, n_atoms
            do j = i + 1, n_atoms
               invd = inv_distance(i, j)
               rij = atom_coords(:, i) - atom_coords(:, j)

               mu = (atom_r(i) - atom_r(j))*invd
               nu = mu + shift(i, j)*(1.0_dp - mu*mu)

               if (scheme == PARTITION_BECKE) then
                  s = becke_cutoff(nu)
                  ds_dnu = becke_cutoff_derivative(nu)
               else
                  s = stratmann_cutoff(nu)
                  ds_dnu = stratmann_cutoff_derivative(nu)
               end if

               dnu_dmu = 1.0_dp - 2.0_dp*shift(i, j)*mu

               ! mu depends on R_i and R_j only: through the distances from the
               ! point, and through the internuclear distance in the denominator.
               dmu_i = -unit_vec(:, i)*invd - (atom_r(i) - atom_r(j))*invd**3*rij
               dmu_j = unit_vec(:, j)*invd + (atom_r(i) - atom_r(j))*invd**3*rij

               ds_i = ds_dnu*dnu_dmu*dmu_i
               ds_j = ds_dnu*dnu_dmu*dmu_j

               ! Product rule, value and derivative together. Order matters:
               ! the derivative uses the cell value from before this factor.
               do a = 1, n_atoms
                  dcell(:, a, i) = dcell(:, a, i)*s
                  dcell(:, a, j) = dcell(:, a, j)*(1.0_dp - s)
               end do
               dcell(:, i, i) = dcell(:, i, i) + cell(i)*ds_i
               dcell(:, j, i) = dcell(:, j, i) + cell(i)*ds_j
               dcell(:, i, j) = dcell(:, i, j) - cell(j)*ds_i
               dcell(:, j, j) = dcell(:, j, j) - cell(j)*ds_j

               cell(i) = cell(i)*s
               cell(j) = cell(j)*(1.0_dp - s)
            end do
         end do

         total = sum(cell)
         if (total <= 0.0_dp) cycle

         do a = 1, n_atoms
            dtotal = 0.0_dp
            do i = 1, n_atoms
               dtotal = dtotal + dcell(:, a, i)
            end do
            dweights(:, a, k) = (dcell(:, a, owner(k))*total &
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
