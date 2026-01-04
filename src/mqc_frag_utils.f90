!! Fragment generation and manipulation utilities
module mqc_frag_utils
   !! Provides combinatorial functions and algorithms for generating molecular
   !! fragments, managing fragment lists, and performing many-body expansion calculations.
   !!
   !! This module re-exports functionality from specialized modules:
   !! - mqc_combinatorics: Pure combinatorial mathematics
   !! - mqc_fragment_lookup: Hash-based fragment index lookup
   !! - mqc_gmbe_utils: GMBE intersection and PIE enumeration
   use pic_types, only: int32, int64, dp, int_index
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_combinatorics, only: &
      binomial, &
      get_nfrags, &
      create_monomer_list, &
      generate_fragment_list, &
      combine, &
      get_next_combination, &
      next_combination_init, &
      next_combination, &
      print_combos, &
      calculate_fragment_distances
   use mqc_fragment_lookup, only: fragment_lookup_t
   use mqc_gmbe_utils, only: &
      find_fragment_intersection, &
      generate_intersections, &
      compute_polymer_atoms, &
      generate_polymer_intersections, &
      gmbe_enumerate_pie_terms
   implicit none
   private

   ! Re-export from mqc_combinatorics
   public :: binomial
   public :: create_monomer_list
   public :: generate_fragment_list
   public :: combine
   public :: get_nfrags
   public :: get_next_combination
   public :: next_combination_init
   public :: next_combination
   public :: print_combos
   public :: calculate_fragment_distances

   ! Re-export from mqc_fragment_lookup
   public :: fragment_lookup_t

   ! Re-export from mqc_gmbe_utils
   public :: find_fragment_intersection
   public :: generate_intersections
   public :: compute_polymer_atoms
   public :: generate_polymer_intersections
   public :: gmbe_enumerate_pie_terms

   ! Local utilities
   public :: apply_distance_screening
   public :: sort_fragments_by_size

contains

   subroutine apply_distance_screening(polymers, total_fragments, sys_geom, driver_config, max_level)
      !! Apply distance-based screening to filter out fragments that exceed cutoff distances
      !! Modifies polymers array in-place and updates total_fragments count
      use mqc_physical_fragment, only: calculate_monomer_distance
      use mqc_config_adapter, only: driver_config_t

      integer, intent(inout) :: polymers(:, :)
      integer(int64), intent(inout) :: total_fragments
      type(system_geometry_t), intent(in) :: sys_geom
      type(driver_config_t), intent(in) :: driver_config
      integer, intent(in) :: max_level

      integer(int64) :: i, fragments_kept
      integer :: fragment_size
      real(dp) :: fragment_distance, cutoff
      integer(int64) :: fragments_screened
      integer :: nmer_level

      ! Check if we have cutoffs to apply
      if (.not. allocated(driver_config%fragment_cutoffs)) then
         return  ! No screening needed
      end if

      fragments_kept = 0_int64
      fragments_screened = 0_int64

      ! Loop through all fragments and filter based on distance
      do i = 1_int64, total_fragments
         fragment_size = count(polymers(i, :) > 0)

         ! Monomers are always kept (distance = 0)
         if (fragment_size == 1) then
            fragments_kept = fragments_kept + 1_int64
            if (fragments_kept /= i) then
               ! Compact array - move this fragment to the kept position
               polymers(fragments_kept, :) = polymers(i, :)
            end if
            cycle
         end if

         ! For n-mers (n >= 2), check if we have a cutoff and apply screening
         nmer_level = fragment_size
         if (nmer_level > size(driver_config%fragment_cutoffs)) then
            ! No cutoff defined for this level - keep the fragment
            fragments_kept = fragments_kept + 1_int64
            if (fragments_kept /= i) then
               polymers(fragments_kept, :) = polymers(i, :)
            end if
            cycle
         end if

         cutoff = driver_config%fragment_cutoffs(nmer_level)
         if (cutoff <= 0.0_dp) then
            ! Negative or zero cutoff means no screening for this level
            fragments_kept = fragments_kept + 1_int64
            if (fragments_kept /= i) then
               polymers(fragments_kept, :) = polymers(i, :)
            end if
            cycle
         end if

         ! Calculate distance for this fragment
         fragment_distance = calculate_monomer_distance(sys_geom, polymers(i, 1:fragment_size))

         ! Apply cutoff filter
         if (fragment_distance <= cutoff) then
            ! Keep this fragment
            fragments_kept = fragments_kept + 1_int64
            if (fragments_kept /= i) then
               polymers(fragments_kept, :) = polymers(i, :)
            end if
         else
            ! Screen out this fragment
            fragments_screened = fragments_screened + 1_int64
         end if
      end do

      ! Update total fragment count
      if (fragments_screened > 0) then
         call logger%info("Distance-based screening applied:")
         call logger%info("  Fragments before screening: "//to_char(total_fragments))
         call logger%info("  Fragments screened out: "//to_char(fragments_screened))
         call logger%info("  Fragments kept: "//to_char(fragments_kept))
         total_fragments = fragments_kept
      end if

   end subroutine apply_distance_screening

   subroutine sort_fragments_by_size(polymers, total_fragments, max_level)
      !! Sort fragments by size (largest first) for better load balancing
      !! Uses in-place sorting to reorder the polymers array
      !! Larger fragments (e.g., tetramers) are computed before smaller ones (e.g., dimers)
      use pic_sorting, only: sort_index

      integer, intent(inout) :: polymers(:, :)
      integer(int64), intent(in) :: total_fragments
      integer, intent(in) :: max_level

      integer(int64), allocatable :: fragment_sizes(:)
      integer(int_index), allocatable :: sort_indices(:)
      integer, allocatable :: polymers_copy(:, :)
      integer(int64) :: i, j, sorted_idx
      integer :: fragment_size

      ! Nothing to sort if we have 1 or fewer fragments
      if (total_fragments <= 1) return

      ! Allocate arrays for sorting (0-indexed for PIC library)
      allocate (fragment_sizes(0:total_fragments - 1))
      allocate (sort_indices(0:total_fragments - 1))

      ! Calculate fragment sizes
      do i = 0, total_fragments - 1
         fragment_size = count(polymers(i + 1, :) > 0)
         fragment_sizes(i) = int(fragment_size, int64)
      end do

      ! Get sort permutation in descending order (largest first)
      call sort_index(fragment_sizes, sort_indices, reverse=.true.)

      ! Reorder polymers array based on sort permutation
      allocate (polymers_copy(size(polymers, 1), size(polymers, 2)))
      polymers_copy = polymers

      ! Reorder: new position j gets data from original position sort_indices(j)
      ! NOTE: sort_indices already contains 1-indexed values, so don't add 1!
      do j = 0, total_fragments - 1
         sorted_idx = sort_indices(j)  ! Already 1-indexed!
         polymers(j + 1, :) = polymers_copy(sorted_idx, :)
      end do

      deallocate (polymers_copy)
      deallocate (fragment_sizes)
      deallocate (sort_indices)

      call logger%info("Fragments sorted by size (largest first) for load balancing")

   end subroutine sort_fragments_by_size

end module mqc_frag_utils
