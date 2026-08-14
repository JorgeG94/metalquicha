module mqc_dft_radial
   !! Radial quadrature for DFT molecular grids
   !!
   !! The Treutler-Ahlrichs M4 mapping [JCP 102, 346 (1995)] takes the Chebyshev
   !! nodes of the second kind on (-1, 1) to (0, infinity):
   !!
   !!    r(x) = -(xi/ln 2) (1 + x)^alpha ln((1 - x)/2)
   !!
   !! with alpha = 0.6 and xi an element-specific scale. The weight returned is
   !! dr/di -- the mapping Jacobian times the Chebyshev weight -- and does *not*
   !! include the r^2 of the volume element or the 4*pi of the sphere.
   !!
   !! That split is deliberate and worth stating, because getting it wrong gives
   !! a grid that integrates to a plausible but wrong number rather than failing.
   !! A caller assembling a molecular grid wants
   !!
   !!    volume weight = 4*pi * r^2 * dr * (angular weight)
   !!
   !! with the angular weight from `mqc_lebedev`, whose weights sum to 1. This is
   !! the same convention PySCF uses, which is what makes the two directly
   !! comparable. `radial_volume_weights` does the r^2 and 4*pi for callers that
   !! want a ready-made spherical quadrature.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_physical_constants, only: PI
   use mqc_dft_radial_data, only: TREUTLER_XI_TABLE, BRAGG_RADII_TABLE, &
                                  MAX_XI_Z, MAX_BRAGG_Z
   implicit none
   private

   public :: treutler_xi              !! Element radius parameter of the M4 mapping
   public :: bragg_radius             !! Bragg-Slater radius, Bohr
   public :: treutler_ahlrichs_radial  !! M4 radial nodes and mapping weights
   public :: radial_volume_weights    !! 4*pi*r^2*dr, for a full spherical quadrature

   !> Exponent of the M4 mapping
   real(dp), parameter :: M4_ALPHA = 0.6_dp

   !> Used where an element has no tabulated value. Every element up to Z=103
   !> has one, so this is only reached by superheavy placeholders.
   real(dp), parameter :: DEFAULT_XI = 1.0_dp
   real(dp), parameter :: DEFAULT_BRAGG = 1.0_dp

contains

   pure function treutler_xi(atomic_number) result(xi)
      !! Element radius parameter of the Treutler-Ahlrichs mapping
      !!
      !! H to Kr are the values optimised in the original paper; heavier
      !! elements come from Psi4. Z=0 is a ghost atom.
      integer, intent(in) :: atomic_number
      real(dp) :: xi

      if (atomic_number >= 0 .and. atomic_number <= MAX_XI_Z) then
         xi = TREUTLER_XI_TABLE(atomic_number)
      else
         xi = DEFAULT_XI
      end if
   end function treutler_xi

   pure function bragg_radius(atomic_number) result(radius)
      !! Bragg-Slater atomic radius in Bohr
      integer, intent(in) :: atomic_number
      real(dp) :: radius

      if (atomic_number >= 0 .and. atomic_number <= MAX_BRAGG_Z) then
         radius = BRAGG_RADII_TABLE(atomic_number)
      else
         radius = DEFAULT_BRAGG
      end if
   end function bragg_radius

   subroutine treutler_ahlrichs_radial(n_points, atomic_number, r, dr, error)
      !! Treutler-Ahlrichs M4 radial mesh for one element
      !!
      !! Nodes come out ascending. `dr` is the mapping weight alone: multiply by
      !! 4*pi*r^2 for a spherical volume element, or use
      !! `radial_volume_weights`.
      integer, intent(in) :: n_points               !! Radial shells
      integer, intent(in) :: atomic_number          !! Z, for the xi parameter
      real(dp), allocatable, intent(out) :: r(:)    !! Nodes, Bohr, ascending
      real(dp), allocatable, intent(out) :: dr(:)   !! Mapping weights
      type(error_t), intent(inout) :: error

      real(dp) :: xi, step, scale, x, sin_x, log_term, map_term
      integer :: i, j

      if (n_points < 1) then
         call error%set(ERROR_VALIDATION, "radial quadrature needs at least one point")
         return
      end if

      allocate (r(n_points), dr(n_points))

      xi = treutler_xi(atomic_number)
      step = PI/real(n_points + 1, dp)
      scale = xi/log(2.0_dp)

      do i = 1, n_points
         x = cos(real(i, dp)*step)
         sin_x = sin(real(i, dp)*step)
         log_term = log((1.0_dp - x)/2.0_dp)
         map_term = (1.0_dp + x)**M4_ALPHA

         ! Filled back to front: i ascending gives x descending, so the nodes
         ! come out descending and are reversed here into ascending order.
         j = n_points - i + 1
         r(j) = -scale*map_term*log_term
         dr(j) = step*sin_x*scale*map_term &
                 *(-M4_ALPHA/(1.0_dp + x)*log_term + 1.0_dp/(1.0_dp - x))
      end do
   end subroutine treutler_ahlrichs_radial

   pure subroutine radial_volume_weights(r, dr, weights)
      !! Turn mapping weights into spherical volume weights, 4*pi*r^2*dr
      !!
      !! Use this only when integrating a spherically symmetric function. For a
      !! molecular grid the 4*pi belongs with the angular weights, which sum to
      !! 1, so applying it here as well would count it twice.
      real(dp), intent(in) :: r(:)         !! Radial nodes
      real(dp), intent(in) :: dr(:)        !! Mapping weights
      real(dp), intent(out) :: weights(:)  !! 4*pi*r^2*dr

      weights = 4.0_dp*PI*r*r*dr
   end subroutine radial_volume_weights

end module mqc_dft_radial
