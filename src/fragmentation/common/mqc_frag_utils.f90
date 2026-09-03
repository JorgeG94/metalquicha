!! Fragment generation and manipulation utilities
module mqc_frag_utils
   !! The term list an expansion runs over: how it is built, screened and
   !! ordered.
   !!
   !! Re-exports three specialised modules under one name:
   !!
   !! - `mqc_combinatorics`: pure combinatorial mathematics
   !! - `mqc_fragment_lookup`: hash-based fragment index lookup
   !! - `mqc_gmbe_utils`: GMBE intersection and PIE enumeration
   use pic_types, only: int32, int64, dp, int_index
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_physical_fragment, only: system_geometry_t
   ! TODO(mqc): `get_next_combination` is named twice in this one `only` list,
   ! on the first continuation line and again eight lines down.
   use mqc_combinatorics, only: fragment_size_of, vmfc_subset_key, real_count_of, get_next_combination, &
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
   public :: generate_mbe_term_list
   public :: sort_fragments_by_size

contains

   subroutine generate_mbe_term_list(sys_geom, driver_config, max_level, polymers, total_fragments)
      !! The term list a plain MBE run evaluates, at this geometry
      !!
      !! Monomers first, then every n-mer up to `max_level`, then distance
      !! screening, the counterpoise rows if any, and the size sort. This is
      !! the same list the driver evaluates, for a caller that needs it in
      !! advance -- an optimization freezing the term list, say.
      !!
      !! Monomers are in the list here, unlike `fraglist_t`, which starts at
      !! pairs: `supplied_terms` is fed straight into the expansion and the
      !! expansion expects them.
      !!
      !! The result is closed under subsets. `fragment_should_be_screened`
      !! drops an n-mer if *any* of its k-subsets exceeds the k-mer cutoff, so
      !! a surviving trimer's dimers all survived too.
      use mqc_config_adapter, only: driver_config_t

      type(system_geometry_t), intent(in) :: sys_geom
      type(driver_config_t), intent(in) :: driver_config
      integer, intent(in) :: max_level
      integer, allocatable, intent(out) :: polymers(:, :)
      integer(int64), intent(out) :: total_fragments

      integer, allocatable :: monomers(:)
      integer(int64) :: n_rows
      integer :: imon

      n_rows = get_nfrags(sys_geom%n_monomers, max_level)

      allocate (monomers(sys_geom%n_monomers))
      allocate (polymers(n_rows, max_level))
      polymers = 0

      call create_monomer_list(monomers)

      total_fragments = 0_int64
      do imon = 1, sys_geom%n_monomers
         total_fragments = total_fragments + 1_int64
         polymers(total_fragments, 1) = imon
      end do

      call generate_fragment_list(monomers, max_level, polymers, total_fragments)
      deallocate (monomers)

      call apply_distance_screening(polymers, total_fragments, sys_geom, driver_config, max_level)

      ! Counterpoise rows are added after screening, so a pair that was screened
      ! out does not bring its ghosted monomers along with it, and before the
      ! size sort, so they are ordered with everything else.
      if (driver_config%counterpoise == "vmfc") then
         call add_vmfc_rows(polymers, total_fragments, max_level)
      end if

      call sort_fragments_by_size(polymers, total_fragments, max_level)

   end subroutine generate_mbe_term_list

   subroutine add_vmfc_rows(polymers, total_fragments, max_level)
      !! Add, for every n-mer, the subfragments solved in that n-mer's basis
      !!
      !! For the pair `[i,j]` that is `[i,-j]` and `[-i,j]`: monomer i with j
      !! present as ghost centres, and the reverse. The pair's correction
      !! subtracts those instead of the bare monomers, so the superposition
      !! error that inflates the pair stands on both sides and cancels.
      !!
      !! The rows are auxiliary -- `is_auxiliary_row` -- and are never summed
      !! into the total.
      integer, allocatable, intent(inout) :: polymers(:, :)
      integer(int64), intent(inout) :: total_fragments
      integer, intent(in) :: max_level

      integer, allocatable :: grown(:, :)
      integer :: parent(max_level), key(max_level), chosen(max_level)
      integer(int64) :: f, n_added, capacity
      integer :: n, k, i
      logical :: has_next

      if (max_level < 2) return

      ! Every n-mer of size n contributes 2^n - 2 proper subsets, so the worst
      ! case is bounded but not small.
      capacity = total_fragments*(2_int64**max_level)
      allocate (grown(capacity, max_level))
      grown = 0
      grown(1:total_fragments, :) = polymers(1:total_fragments, :)
      n_added = total_fragments

      do f = 1, total_fragments
         n = int(real_count_of(polymers(f, :)))
         if (n < 2) cycle              ! a monomer has no subsets to ghost
         parent(1:n) = polymers(f, 1:n)

         do k = 1, n - 1
            do i = 1, k
               chosen(i) = i
            end do
            do
               call vmfc_subset_key(parent(1:n), n, chosen(1:k), k, key(1:n))
               n_added = n_added + 1_int64
               grown(n_added, 1:n) = key(1:n)
               call get_next_combination(chosen, k, n, has_next)
               if (.not. has_next) exit
            end do
         end do
      end do

      call move_alloc(grown, polymers)
      total_fragments = n_added
   end subroutine add_vmfc_rows

   subroutine apply_distance_screening(polymers, total_fragments, sys_geom, driver_config, max_level)
      !! Drop the fragments beyond their level's cutoff, in place
      !!
      !! `polymers` is compacted and `total_fragments` reduced. An n-mer goes
      !! if *any* of its k-subsets exceeds the k-mer cutoff, so the surviving
      !! list stays closed under subsets -- `compute_mbe` looks those subsets
      !! up and fails on a missing one.
      use mqc_physical_fragment, only: calculate_monomer_distance
      use mqc_config_adapter, only: driver_config_t

      integer, intent(inout) :: polymers(:, :)
      integer(int64), intent(inout) :: total_fragments
      type(system_geometry_t), intent(in) :: sys_geom
      type(driver_config_t), intent(in) :: driver_config
      integer, intent(in) :: max_level

      integer(int64) :: i, fragments_kept
      integer :: fragment_size
      integer(int64) :: fragments_screened
      logical :: should_screen

      ! Check if we have cutoffs to apply
      if (.not. allocated(driver_config%fragment_cutoffs)) then
         return  ! No screening needed
      end if

      fragments_kept = 0_int64
      fragments_screened = 0_int64

      ! Loop through all fragments and filter based on distance
      do i = 1_int64, total_fragments
         fragment_size = fragment_size_of(polymers(i, :))

         ! Monomers are always kept (distance = 0)
         if (fragment_size == 1) then
            fragments_kept = fragments_kept + 1_int64
            if (fragments_kept /= i) then
               ! Compact array - move this fragment to the kept position
               polymers(fragments_kept, :) = polymers(i, :)
            end if
            cycle
         end if

         ! For n-mers (n >= 2), check if this fragment or any of its subsets should be screened
         should_screen = fragment_should_be_screened(polymers(i, 1:fragment_size), fragment_size, &
                                                     sys_geom, driver_config)

         if (.not. should_screen) then
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

   function fragment_should_be_screened(fragment, n, sys_geom, driver_config) result(should_screen)
      !! Whether this fragment or any of its k-subsets, `k >= 2`, exceeds the
      !! k-mer cutoff
      use mqc_physical_fragment, only: calculate_monomer_distance
      use mqc_config_adapter, only: driver_config_t

      integer, intent(in) :: fragment(:)
      integer, intent(in) :: n
      type(system_geometry_t), intent(in) :: sys_geom
      type(driver_config_t), intent(in) :: driver_config
      logical :: should_screen

      integer :: subset_size, num_cutoffs
      integer :: indices(n), subset(n)
      integer :: j
      real(dp) :: distance, cutoff
      logical :: has_next

      should_screen = .false.
      num_cutoffs = size(driver_config%fragment_cutoffs)

      ! Check all subset sizes from 2 up to n (the full fragment)
      ! If any k-subset exceeds the k-mer cutoff, screen this fragment
      do subset_size = 2, n
         ! Skip if no cutoff defined for this level
         if (subset_size > num_cutoffs) cycle

         cutoff = driver_config%fragment_cutoffs(subset_size)
         ! Skip if cutoff is non-positive (no screening for this level)
         if (cutoff <= 0.0_dp) cycle

         ! Initialize first combination indices
         do j = 1, subset_size
            indices(j) = j
         end do

         ! Loop through all combinations of this size
         do
            ! Build current subset
            do j = 1, subset_size
               subset(j) = fragment(indices(j))
            end do

            ! Calculate distance for this subset
            distance = calculate_monomer_distance(sys_geom, subset(1:subset_size))

            ! If subset exceeds cutoff, screen the whole fragment
            if (distance > cutoff) then
               should_screen = .true.
               return
            end if

            ! Get next combination
            call get_next_combination(indices, subset_size, n, has_next)
            if (.not. has_next) exit
         end do
      end do

   end function fragment_should_be_screened

   ! cannot make this pure because sort is not pure
   subroutine sort_fragments_by_size(polymers, total_fragments, max_level)
      !! Reorder `polymers` in place, largest fragment first
      !!
      !! The expensive terms then start first, which is what balances the load
      !! across ranks.
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
         fragment_size = fragment_size_of(polymers(i + 1, :))
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

      call logger%info("Fragments queue sorted!")

   end subroutine sort_fragments_by_size

end module mqc_frag_utils
