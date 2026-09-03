module mqc_dft_prune
   !! Angular pruning: fewer Lebedev points where the density does not need them
   !!
   !! The NWChem scheme divides each atom's radial range into five zones by
   !! r/R_bragg, with thresholds depending on whether the element is H/He,
   !! second row, or heavier, and assigns orders
   !!
   !!    50, 86, one below the target, the target, one below the target
   !!
   !! reading outwards. Pruning is an approximation and changes the answer,
   !! within the accuracy the grid claims.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_dft_radial, only: bragg_radius
   use mqc_dft_partition, only: TINY_RADIUS
   implicit none
   private

   public :: PRUNE_NONE, PRUNE_NWCHEM
   public :: prune_angular_orders
   public :: prune_scheme_name

   integer, parameter :: PRUNE_NONE = 0    !! Same order at every radius
   integer, parameter :: PRUNE_NWCHEM = 1  !! NWChem five-zone scheme

   integer, parameter :: N_ZONES = 5
   integer, parameter :: N_THRESHOLDS = 4

   integer, parameter :: N_PRUNE_ORDERS = 29
   integer, parameter :: PRUNE_ORDERS(N_PRUNE_ORDERS) = &
                         [38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, &
                          974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, &
                          5294, 5810]
      !! Lebedev orders the scheme may select: the sequence from 38 upwards, not
      !! from 6. The scheme indexes into it, so the three smallest grids are
      !! never chosen.

   real(dp), parameter :: ALPHAS(N_THRESHOLDS, 3) = reshape([ &
                                                            0.25_dp, 0.5_dp, 1.0_dp, 4.5_dp, &
                                                            0.1667_dp, 0.5_dp, 0.9_dp, 3.5_dp, &
                                                            0.1_dp, 0.4_dp, 0.8_dp, 2.5_dp], &
                                                            [N_THRESHOLDS, 3])
      !! Zone boundaries in r/R_bragg, one row per element class

   integer, parameter :: MIN_PRUNABLE = 50
      !! Below this order there is nothing to gain, so pruning is skipped

contains

   pure function prune_scheme_name(scheme) result(name)
      !! Human-readable scheme name
      integer, intent(in) :: scheme
      character(len=:), allocatable :: name

      select case (scheme)
      case (PRUNE_NONE)
         name = "none"
      case (PRUNE_NWCHEM)
         name = "nwchem"
      case default
         name = "unknown"
      end select
   end function prune_scheme_name

   subroutine prune_angular_orders(scheme, atomic_number, r, n_angular, orders, error)
      !! Lebedev order for each radial shell of one atom
      integer, intent(in) :: scheme             !! PRUNE_NONE or PRUNE_NWCHEM
      integer, intent(in) :: atomic_number      !! Z, for the radius and element class
      real(dp), intent(in) :: r(:)              !! Radial nodes, Bohr
      integer, intent(in) :: n_angular          !! Target order for the valence zone
      integer, intent(out) :: orders(:)         !! Order per shell, same size as r
      type(error_t), intent(inout) :: error

      integer :: zone_order(N_ZONES)
      integer :: target_index, element_class, i, zone
      real(dp) :: scaled

      if (size(orders) /= size(r)) then
         call error%set(ERROR_VALIDATION, "prune: orders must match the radial mesh")
         return
      end if

      if (scheme == PRUNE_NONE .or. n_angular < MIN_PRUNABLE) then
         orders = n_angular
         return
      end if

      if (scheme /= PRUNE_NWCHEM) then
         call error%set(ERROR_VALIDATION, "prune: unknown scheme")
         return
      end if

      if (n_angular == MIN_PRUNABLE) then
         ! The smallest prunable grid has no room to step down twice, so the
         ! scheme steps up in the middle instead of down at the edges.
         zone_order = [PRUNE_ORDERS(2), PRUNE_ORDERS(3), PRUNE_ORDERS(3), &
                       PRUNE_ORDERS(3), PRUNE_ORDERS(2)]
      else
         target_index = order_index(n_angular)
         if (target_index == 0) then
            call error%set(ERROR_VALIDATION, "prune: target order is not a Lebedev order")
            return
         end if
         zone_order = [PRUNE_ORDERS(2), PRUNE_ORDERS(4), PRUNE_ORDERS(target_index - 1), &
                       PRUNE_ORDERS(target_index), PRUNE_ORDERS(target_index - 1)]
      end if

      if (atomic_number <= 2) then
         element_class = 1
      else if (atomic_number <= 10) then
         element_class = 2
      else
         element_class = 3
      end if

      do i = 1, size(r)
         scaled = r(i)/(bragg_radius(atomic_number) + TINY_RADIUS)
         ! Zone is the number of thresholds the scaled radius has passed.
         zone = count(scaled > ALPHAS(:, element_class)) + 1
         orders(i) = zone_order(zone)
      end do
   end subroutine prune_angular_orders

   pure function order_index(n_angular) result(idx)
      !! Position of an order in the prunable sequence, or 0 if absent
      integer, intent(in) :: n_angular
      integer :: idx
      integer :: i

      idx = 0
      do i = 1, N_PRUNE_ORDERS
         if (PRUNE_ORDERS(i) == n_angular) then
            idx = i
            exit
         end if
      end do
   end function order_index

end module mqc_dft_prune
