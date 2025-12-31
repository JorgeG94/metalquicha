!! Fragment generation and manipulation utilities
module mqc_frag_utils
   !! Provides combinatorial functions and algorithms for generating molecular
   !! fragments, managing fragment lists, and performing many-body expansion calculations.
   use pic_types, only: default_int, int32, int64, dp
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use pic_sorting, only: sort
   use pic_hash_32bit, only: fnv_1a_hash
   implicit none
   private
   public :: binomial              !! Binomial coefficient calculation
   public :: create_monomer_list   !! Generate sequential monomer indices
   public :: generate_fragment_list  !! Generate all fragments up to max level
   public :: combine               !! Generate all combinations of size r (helper for single-level generation)
   public :: get_nfrags            !! Calculate total number of fragments
   public :: get_next_combination      !! Generate next combination in sequence
   public :: fragment_lookup_t     !! Hash-based lookup table for fast fragment index retrieval
   public :: find_fragment_intersection  !! Find shared atoms between two fragments (for GMBE)
   public :: generate_intersections     !! Generate all k-way intersections for GMBE (base fragments)
   public :: compute_polymer_atoms      !! Compute atom list for a polymer (union of base fragments)
   public :: generate_polymer_intersections  !! Generate intersections for polymers at any level
   public :: gmbe_enumerate_pie_terms   !! DFS-based PIE coefficient enumeration for GMBE(N)

   type :: hash_entry_t
      !! Single entry in hash table (private helper type)
      integer, allocatable :: key(:)      !! Sorted monomer indices
      integer(int64) :: value             !! Fragment index
      type(hash_entry_t), pointer :: next => null()  !! Chain for collisions
   end type hash_entry_t

   type :: fragment_lookup_t
      !! Hash-based lookup table for O(1) fragment index retrieval
      integer :: table_size = 0
      type(hash_entry_t), allocatable :: table(:)
      integer(int64) :: n_entries = 0
      logical :: initialized = .false.
   contains
      procedure :: init => fragment_lookup_init
      procedure :: insert => fragment_lookup_insert
      procedure :: find => fragment_lookup_find
      procedure :: destroy => fragment_lookup_destroy
   end type fragment_lookup_t

contains

   pure function get_nfrags(n_monomers, max_level) result(n_expected_fragments)
      !! Calculate total number of fragments for given system size and max level
      !!
      !! Computes the sum of binomial coefficients C(n,k) for k=1 to max_level,
      !! representing all possible fragments from monomers to max_level-mers.
      !! Uses int64 to handle large fragment counts that overflow int32.
      integer(default_int), intent(in) :: n_monomers  !! Number of monomers in system
      integer(default_int), intent(in) :: max_level   !! Maximum fragment size
      integer(int64) :: n_expected_fragments     !! Total fragment count
      integer(default_int) :: i  !! Loop counter

      n_expected_fragments = 0_int64
      do i = 1, max_level
         n_expected_fragments = n_expected_fragments + binomial(n_monomers, i)
      end do
   end function get_nfrags

   pure function binomial(n, r) result(c)
      !! Compute binomial coefficient C(n,r) = n! / (r! * (n-r)!)
      !!
      !! Calculates "n choose r" using iterative algorithm to avoid
      !! factorial overflow for large numbers.
      !! Uses int64 to handle large combinatorial values that overflow int32.
      integer(default_int), intent(in) :: n  !! Total number of items
      integer(default_int), intent(in) :: r  !! Number of items to choose
      integer(int64) :: c              !! Binomial coefficient result
      integer(default_int) :: i              !! Loop counter

      if (r == 0 .or. r == n) then
         c = 1_int64
      else if (r > n) then
         c = 0_int64
      else
         c = 1_int64
         do i = 1, r
            c = c*int(n - i + 1, int64)/int(i, int64)
         end do
      end if
   end function binomial

   pure subroutine create_monomer_list(monomers)
      !! Generate a list of monomer indices from 1 to N
      integer(default_int), allocatable, intent(inout) :: monomers(:)
      integer(default_int) :: i, length

      length = size(monomers, 1)

      do i = 1, length
         monomers(i) = i
      end do

   end subroutine create_monomer_list

   recursive subroutine generate_fragment_list(monomers, max_level, polymers, count)
      !! Generate all possible fragments (combinations of monomers) up to max_level
      !! Uses int64 for count to handle large numbers of fragments that overflow int32.
      integer(default_int), intent(in) :: monomers(:), max_level
      integer(default_int), intent(inout) :: polymers(:, :)
      integer(int64), intent(inout) :: count
      integer(default_int) :: r, n

      n = size(monomers, 1)
      do r = 2, max_level
         call combine(monomers, n, r, polymers, count)
      end do
   end subroutine generate_fragment_list

   recursive subroutine combine(arr, n, r, out_array, count)
      !! Generate all combinations of size r from array arr of size n
      !! Uses int64 for count to handle large numbers of combinations that overflow int32.
      integer(default_int), intent(in) :: arr(:)
      integer(default_int), intent(in) :: n, r
      integer(default_int), intent(inout) :: out_array(:, :)
      integer(int64), intent(inout) :: count
      integer(default_int) :: data(r)
      call combine_util(arr, n, r, 1, data, 1, out_array, count)
   end subroutine combine

   recursive subroutine combine_util(arr, n, r, index, data, i, out_array, count)
      !! Utility for generating combinations recursively
      !! Uses int64 for count to handle large numbers of combinations that overflow int32.
      integer(default_int), intent(in) :: arr(:), n, r, index, i
      integer(default_int), intent(inout) :: data(:), out_array(:, :)
      integer(int64), intent(inout) :: count
      integer(default_int) :: j

      if (index > r) then
         count = count + 1_int64
         out_array(count, 1:r) = data(1:r)
         return
      end if

      do j = i, n
         data(index) = arr(j)
         call combine_util(arr, n, r, index + 1, data, j + 1, out_array, count)
      end do
   end subroutine combine_util

   subroutine print_combos(out_array, count, max_len)
      !! Print combinations stored in out_array
      !! Uses int64 for count to handle large numbers of combinations that overflow int32.
      integer(default_int), intent(in) :: out_array(:, :), max_len
      integer(int64), intent(in) :: count
      integer(int64) :: i
      integer(default_int) :: j

      do i = 1_int64, count
         do j = 1, max_len
            if (out_array(i, j) == 0) exit
            write (*, '(I0)', advance='no') out_array(i, j)
            if (j < max_len .and. out_array(i, j + 1) /= 0) then
               write (*, '(A)', advance='no') ":"
            end if
         end do
         write (*, *)  ! newline
      end do
   end subroutine print_combos

   pure subroutine get_next_combination(indices, k, n, has_next)
   !! Generate next combination (updates indices in place)
   !! has_next = .true. if there's a next combination
      integer, intent(inout) :: indices(:)
      integer, intent(in)    :: k, n
      logical, intent(out)   :: has_next
      integer                :: i

      has_next = .true.

      i = k
      do while (i >= 1)
         if (indices(i) < n - k + i) then
            indices(i) = indices(i) + 1
            do while (i < k)
               i = i + 1
               indices(i) = indices(i - 1) + 1
            end do
            return
         end if
         i = i - 1
      end do

      has_next = .false.
   end subroutine get_next_combination

   subroutine fragment_lookup_init(this, estimated_entries)
      !! Initialize hash table with estimated size
      class(fragment_lookup_t), intent(inout) :: this
      integer(int64), intent(in) :: estimated_entries

      integer :: i

      ! Use prime number close to estimated size for better distribution
      this%table_size = next_prime_internal(int(estimated_entries*1.3_dp))
      allocate (this%table(this%table_size))

      ! Initialize all entries as empty
      do i = 1, this%table_size
         nullify (this%table(i)%next)
      end do

      this%n_entries = 0
      this%initialized = .true.
   end subroutine fragment_lookup_init

   subroutine fragment_lookup_insert(this, monomers, n, fragment_idx)
      !! Insert a monomer combination -> fragment index mapping
      class(fragment_lookup_t), intent(inout) :: this
      integer, intent(in) :: monomers(:), n
      integer(int64), intent(in) :: fragment_idx

      integer(int32) :: hash_val
      integer :: bucket
      type(hash_entry_t), pointer :: new_entry
      integer, allocatable :: sorted_key(:)

      if (.not. this%initialized) error stop "Hash table not initialized"

      ! Sort monomers for canonical key
      allocate (sorted_key(n))
      sorted_key = monomers(1:n)
      call sort(sorted_key)

      ! Compute hash
      hash_val = fnv_1a_hash(sorted_key)
      bucket = 1 + modulo(hash_val, int(this%table_size, int32))

      ! Check if this is the first entry in bucket
      if (.not. allocated(this%table(bucket)%key)) then
         ! First entry in this bucket - use the head entry
         allocate (this%table(bucket)%key(n))
         this%table(bucket)%key = sorted_key
         this%table(bucket)%value = fragment_idx
         this%n_entries = this%n_entries + 1
      else
         ! Bucket already has entries - chain new entry
         allocate (new_entry)
         allocate (new_entry%key(n))
         new_entry%key = sorted_key
         new_entry%value = fragment_idx
         new_entry%next => this%table(bucket)%next
         this%table(bucket)%next => new_entry
         this%n_entries = this%n_entries + 1
      end if

      deallocate (sorted_key)
   end subroutine fragment_lookup_insert

   function fragment_lookup_find(this, monomers, n) result(idx)
      !! Find fragment index for given monomer combination
      class(fragment_lookup_t), intent(in) :: this
      integer, intent(in) :: monomers(:), n
      integer(int64) :: idx

      integer(int32) :: hash_val
      integer :: bucket, sorted_key(n)
      type(hash_entry_t), pointer :: entry

      ! Sort monomers for canonical key
      sorted_key = monomers(1:n)
      call sort(sorted_key)

      ! Compute hash
      hash_val = fnv_1a_hash(sorted_key)
      bucket = 1 + modulo(hash_val, int(this%table_size, int32))

      ! Search chain
      if (allocated(this%table(bucket)%key)) then
         if (arrays_equal_internal(this%table(bucket)%key, sorted_key, n)) then
            idx = this%table(bucket)%value
            return
         end if
         entry => this%table(bucket)%next
         do while (associated(entry))
            if (arrays_equal_internal(entry%key, sorted_key, n)) then
               idx = entry%value
               return
            end if
            entry => entry%next
         end do
      end if

      ! Not found
      idx = -1
   end function fragment_lookup_find

   subroutine fragment_lookup_destroy(this)
      !! Clean up hash table and all chains
      class(fragment_lookup_t), intent(inout) :: this
      integer :: i
      type(hash_entry_t), pointer :: entry, next_entry

      if (.not. this%initialized) return

      do i = 1, this%table_size
         ! Free chain
         entry => this%table(i)%next
         do while (associated(entry))
            next_entry => entry%next
            if (allocated(entry%key)) deallocate (entry%key)
            deallocate (entry)
            entry => next_entry
         end do
         ! Free bucket head
         if (allocated(this%table(i)%key)) deallocate (this%table(i)%key)
      end do

      deallocate (this%table)
      this%initialized = .false.
   end subroutine fragment_lookup_destroy

   ! Helper functions for hash table
   pure function arrays_equal_internal(a, b, n) result(equal)
      !! Check if two arrays are equal
      integer, intent(in) :: a(:), b(:), n
      logical :: equal
      integer :: i

      equal = .true.
      if (size(a) /= n .or. size(b) /= n) then
         equal = .false.
         return
      end if

      do i = 1, n
         if (a(i) /= b(i)) then
            equal = .false.
            return
         end if
      end do
   end function arrays_equal_internal

   pure function next_prime_internal(n) result(p)
      !! Find next prime number >= n (simple implementation)
      integer, intent(in) :: n
      integer :: p, i
      logical :: is_prime

      p = max(n, 2)
      if (modulo(p, 2) == 0) p = p + 1

      do
         is_prime = .true.
         do i = 3, int(sqrt(real(p))) + 1, 2
            if (modulo(p, i) == 0) then
               is_prime = .false.
               exit
            end if
         end do
         if (is_prime) return
         p = p + 2
      end do
   end function next_prime_internal

   function find_fragment_intersection(frag1_atoms, n1, frag2_atoms, n2, &
                                       intersection, n_intersect) result(has_intersection)
      !! Find shared atoms between two fragments (for GMBE with overlapping fragments)
      !!
      !! This function identifies atoms that appear in both fragments, which is essential
      !! for computing intersection-corrected energies in GMBE.
      !!
      !! Algorithm: O(n1 * n2) brute-force comparison
      !! - Loop through all atoms in fragment 1
      !! - For each atom, check if it appears in fragment 2
      !! - Collect all shared atoms
      !!
      !! Returns:
      !!   .true. if fragments share at least one atom, .false. otherwise
      !!
      !! Output:
      !!   intersection - allocatable array containing shared atom indices
      !!   n_intersect - number of shared atoms
      integer, intent(in) :: frag1_atoms(:)  !! Atom indices in fragment 1 (0-indexed)
      integer, intent(in) :: n1              !! Number of atoms in fragment 1
      integer, intent(in) :: frag2_atoms(:)  !! Atom indices in fragment 2 (0-indexed)
      integer, intent(in) :: n2              !! Number of atoms in fragment 2
      integer, allocatable, intent(out) :: intersection(:)  !! Shared atom indices
      integer, intent(out) :: n_intersect    !! Number of shared atoms
      logical :: has_intersection

      integer :: i, j
      integer, allocatable :: temp_intersection(:)
      integer :: temp_count

      ! Allocate temporary array (max possible size is min(n1, n2))
      allocate (temp_intersection(min(n1, n2)))
      temp_count = 0

      ! Find all shared atoms
      do i = 1, n1
         do j = 1, n2
            if (frag1_atoms(i) == frag2_atoms(j)) then
               ! Found a shared atom
               temp_count = temp_count + 1
               temp_intersection(temp_count) = frag1_atoms(i)
               exit  ! Move to next atom in frag1
            end if
         end do
      end do

      ! Set output
      n_intersect = temp_count
      has_intersection = (temp_count > 0)

      ! Allocate and copy result if intersection exists
      if (has_intersection) then
         allocate (intersection(n_intersect))
         intersection = temp_intersection(1:n_intersect)
      end if

      deallocate (temp_intersection)

   end function find_fragment_intersection

   subroutine generate_intersections(sys_geom, monomers, polymers, n_monomers, max_intersection_level, &
                                     intersections, intersection_sets, intersection_levels, n_intersections)
      !! Generate all k-way intersections for k=2 to min(max_intersection_level, n_monomers)
      !!\n!! For a system with overlapping fragments, this computes k-way intersections
      !! following the inclusion-exclusion principle for GMBE.
      !! The max_intersection_level parameter controls the maximum depth to avoid combinatorial explosion.
      !!\n!! Algorithm:
      !! - For each k from 2 to min(max_intersection_level, n_monomers):
      !!   - Generate all C(n_monomers, k) combinations
      !!   - For each combination, compute intersection of all k fragments
      !!   - Store non-empty intersections with their level k
      use mqc_physical_fragment, only: system_geometry_t
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomers(:)           !! Monomer indices
      integer, intent(inout) :: polymers(:, :)     !! Output: monomers stored here
      integer, intent(in) :: n_monomers            !! Number of monomers
      integer, intent(in) :: max_intersection_level  !! Maximum k-way intersection depth
      integer, allocatable, intent(out) :: intersections(:, :)  !! Intersection atom lists
      integer, allocatable, intent(out) :: intersection_sets(:, :)  !! Which k-tuple created each intersection
      integer, allocatable, intent(out) :: intersection_levels(:)  !! Level (k) of each intersection
      integer, intent(out) :: n_intersections      !! Number of intersections found

      ! Temporaries for storing intersections
      integer, allocatable :: temp_intersections(:, :)
      integer, allocatable :: temp_sets(:, :)
      integer, allocatable :: temp_levels(:)
      integer, allocatable :: temp_intersection(:)
      integer, allocatable :: current_intersection(:)
      integer :: temp_n_intersect, current_n_intersect
      logical :: has_intersection
      integer :: k, intersection_count, max_atoms, max_intersections, max_k_level
      integer :: i, idx
      integer, allocatable :: combination(:)

      ! Store monomers in polymers array
      polymers(1:n_monomers, 1) = monomers(1:n_monomers)

      if (n_monomers < 2) then
         n_intersections = 0
         return
      end if

      ! Count maximum possible intersections: sum of C(n,k) for k=2 to n
      ! For small n, this is 2^n - n - 1
      max_intersections = 2**n_monomers - n_monomers - 1

      ! Find maximum atoms in any fragment for allocation
      max_atoms = maxval(sys_geom%fragment_sizes(1:n_monomers))

      ! Allocate temporary arrays
      allocate (temp_intersections(max_atoms, max_intersections))
      allocate (temp_sets(n_monomers, max_intersections))
      allocate (temp_levels(max_intersections))
      temp_intersections = 0
      temp_sets = 0
      intersection_count = 0

      ! Determine actual maximum intersection level to use
      max_k_level = min(max_intersection_level, n_monomers)

      if (max_k_level < n_monomers) then
         call logger%info("Generating k-way intersections up to k="//to_char(max_k_level)// &
                          " (limited by max_intersection_level)")
      else
         call logger%info("Generating all k-way intersections for GMBE (inclusion-exclusion principle)")
      end if

      ! Loop over intersection levels k from 2 to max_k_level
      do k = 2, max_k_level
         ! Generate all C(n_monomers, k) combinations
         allocate (combination(k))
         call generate_k_way_intersections_for_level(sys_geom, monomers, n_monomers, k, &
                                                     combination, max_atoms, &
                                                     temp_intersections, temp_sets, temp_levels, intersection_count)
         deallocate (combination)
      end do

      n_intersections = intersection_count

      ! Allocate output arrays
      if (n_intersections > 0) then
         allocate (intersections(max_atoms, n_intersections))
         allocate (intersection_sets(n_monomers, n_intersections))
         allocate (intersection_levels(n_intersections))
         intersections = temp_intersections(1:max_atoms, 1:n_intersections)
         intersection_sets = temp_sets(1:n_monomers, 1:n_intersections)
         intersection_levels = temp_levels(1:n_intersections)

         call logger%info("Generated "//to_char(n_intersections)//" total intersections:")
         do k = 2, max_k_level
            idx = count(intersection_levels == k)
            if (idx > 0) then
               call logger%info("  "//to_char(idx)//" intersections at level "//to_char(k))
            end if
         end do
      else
         call logger%info("No intersections found (fragments are non-overlapping)")
      end if

      deallocate (temp_intersections, temp_sets, temp_levels)

   end subroutine generate_intersections

   recursive subroutine generate_k_way_intersections_for_level(sys_geom, monomers, n_monomers, k, &
                                                               combination, max_atoms, &
                                                           temp_intersections, temp_sets, temp_levels, intersection_count)
      !! Helper to generate all k-way intersections at a specific level k
      use mqc_physical_fragment, only: system_geometry_t
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomers(:), n_monomers, k, max_atoms
      integer, intent(inout) :: combination(:)
      integer, intent(inout) :: temp_intersections(:, :), temp_sets(:, :), temp_levels(:)
      integer, intent(inout) :: intersection_count

      integer, allocatable :: current_intersection(:), temp_intersection(:)
      integer :: current_n_intersect, temp_n_intersect
      integer :: i, j, mono_idx, frag_size
      logical :: has_intersection

      ! Generate combinations using an iterative approach
      call next_combination_init(combination, k)

      do
         ! Compute intersection of all fragments in this combination
         has_intersection = .false.

         ! Start with the first fragment in the combination
         mono_idx = monomers(combination(1))
         frag_size = sys_geom%fragment_sizes(mono_idx)
         allocate (current_intersection(frag_size))
         current_intersection = sys_geom%fragment_atoms(1:frag_size, mono_idx)
         current_n_intersect = frag_size

         ! Intersect with each subsequent fragment
         do i = 2, k
            mono_idx = monomers(combination(i))
            frag_size = sys_geom%fragment_sizes(mono_idx)

            has_intersection = find_fragment_intersection( &
                               current_intersection, current_n_intersect, &
                               sys_geom%fragment_atoms(1:frag_size, mono_idx), frag_size, &
                               temp_intersection, temp_n_intersect)

            if (.not. has_intersection) then
               ! Intersection is empty, break early
               deallocate (current_intersection)
               if (allocated(temp_intersection)) deallocate (temp_intersection)
               exit
            end if

            ! Replace current_intersection with the new intersection
            deallocate (current_intersection)
            allocate (current_intersection(temp_n_intersect))
            current_intersection = temp_intersection
            current_n_intersect = temp_n_intersect
            deallocate (temp_intersection)
         end do

         ! If we have a non-empty intersection, store it
         if (has_intersection .and. allocated(current_intersection)) then
            intersection_count = intersection_count + 1
            temp_intersections(1:current_n_intersect, intersection_count) = current_intersection
            temp_sets(1:k, intersection_count) = monomers(combination)
            temp_levels(intersection_count) = k
            deallocate (current_intersection)
         end if

         ! Get next combination
         if (.not. next_combination(combination, k, n_monomers)) exit
      end do

   end subroutine generate_k_way_intersections_for_level

   subroutine next_combination_init(combination, k)
      !! Initialize combination to [1, 2, ..., k]
      integer, intent(inout) :: combination(:)
      integer, intent(in) :: k
      integer :: i
      do i = 1, k
         combination(i) = i
      end do
   end subroutine next_combination_init

   function next_combination(combination, k, n) result(has_next)
      !! Generate next combination in lexicographic order
      !! Returns .true. if there's a next combination, .false. if we've exhausted all
      integer, intent(inout) :: combination(:)
      integer, intent(in) :: k, n
      logical :: has_next
      integer :: i

      has_next = .true.

      ! Find the rightmost element that can be incremented
      i = k
      do while (i >= 1)
         if (combination(i) < n - k + i) then
            combination(i) = combination(i) + 1
            ! Reset all elements to the right
            do while (i < k)
               i = i + 1
               combination(i) = combination(i - 1) + 1
            end do
            return
         end if
         i = i - 1
      end do

      ! No more combinations
      has_next = .false.

   end function next_combination

   subroutine compute_polymer_atoms(sys_geom, polymer, polymer_size, atom_list, n_atoms)
      !! Compute the atom list for a polymer (union of atoms from base fragments)
      !! polymer(:) contains base fragment indices (1-based)
      use mqc_physical_fragment, only: system_geometry_t
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: polymer(:)  !! Base fragment indices in this polymer
      integer, intent(in) :: polymer_size  !! Number of base fragments in polymer
      integer, allocatable, intent(out) :: atom_list(:)  !! Unique atoms in this polymer
      integer, intent(out) :: n_atoms  !! Number of unique atoms

      integer, allocatable :: temp_atoms(:)
      integer :: i, j, frag_idx, frag_size, temp_count
      logical :: already_present

      ! Allocate temporary array (worst case: all atoms from all fragments)
      allocate (temp_atoms(sys_geom%total_atoms))
      temp_count = 0

      ! Loop through each base fragment in the polymer
      do i = 1, polymer_size
         frag_idx = polymer(i)
         if (frag_idx == 0) exit  ! Padding zeros

         frag_size = sys_geom%fragment_sizes(frag_idx)

         ! Add each atom from this fragment (avoid duplicates)
         do j = 1, frag_size
            already_present = .false.

            ! Check if this atom is already in our list
            block
               integer :: k, current_atom
               current_atom = sys_geom%fragment_atoms(j, frag_idx)
               do k = 1, temp_count
                  if (temp_atoms(k) == current_atom) then
                     already_present = .true.
                     exit
                  end if
               end do
            end block

            ! Add if not already present
            if (.not. already_present) then
               temp_count = temp_count + 1
               temp_atoms(temp_count) = sys_geom%fragment_atoms(j, frag_idx)
            end if
         end do
      end do

      ! Copy to output array
      n_atoms = temp_count
      allocate (atom_list(n_atoms))
      atom_list = temp_atoms(1:n_atoms)

      deallocate (temp_atoms)
   end subroutine compute_polymer_atoms

   subroutine generate_polymer_intersections(sys_geom, polymers, n_polymers, max_level, &
                                             intersections, intersection_sets, intersection_levels, n_intersections)
      !! Generate all k-way intersections for polymers at any level (GMBE-N)
      !! This works with dynamically generated polymers, not just base fragments
      use mqc_physical_fragment, only: system_geometry_t
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: polymers(:, :)  !! Polymer definitions (n_polymers, max_level)
      integer, intent(in) :: n_polymers, max_level
      integer, allocatable, intent(out) :: intersections(:, :)
      integer, allocatable, intent(out) :: intersection_sets(:, :)
      integer, allocatable, intent(out) :: intersection_levels(:)
      integer, intent(out) :: n_intersections

      integer, allocatable :: polymer_atoms(:, :)  !! Atom lists for each polymer
      integer, allocatable :: polymer_n_atoms(:)   !! Number of atoms in each polymer
      integer :: max_atoms_per_polymer
      integer :: i, polymer_size, max_intersection_level

      call logger%info("Computing atom compositions for "//to_char(n_polymers)//" polymers...")

      ! First, compute atom list for each polymer
      ! Find max atoms needed
      max_atoms_per_polymer = 0
      do i = 1, n_polymers
         polymer_size = count(polymers(i, :) > 0)
         ! Worst case: all atoms from all fragments in this polymer
         max_atoms_per_polymer = max(max_atoms_per_polymer, polymer_size*maxval(sys_geom%fragment_sizes))
      end do

      allocate (polymer_atoms(max_atoms_per_polymer, n_polymers))
      allocate (polymer_n_atoms(n_polymers))
      polymer_atoms = 0

      ! Compute atoms for each polymer
      do i = 1, n_polymers
         polymer_size = count(polymers(i, :) > 0)
         block
            integer, allocatable :: atom_list(:)
            integer :: n_atoms
            call compute_polymer_atoms(sys_geom, polymers(i, 1:polymer_size), polymer_size, atom_list, n_atoms)
            polymer_n_atoms(i) = n_atoms
            polymer_atoms(1:n_atoms, i) = atom_list
            deallocate (atom_list)
         end block
      end do

      call logger%info("Finding intersections between polymers...")

      ! For GMBE(N), limit intersections to level N+1 to prevent combinatorial explosion
      ! GMBE(2): dimers → 3-way intersections max
      ! GMBE(3): trimers → 4-way intersections max
      max_intersection_level = max_level + 1
      call logger%info("Limiting intersections to level "//to_char(max_intersection_level)// &
                       " (polymer level "//to_char(max_level)//" + 1)")

      ! Now generate intersections between these polymer atom sets
      call generate_intersections_from_atom_lists(polymer_atoms, polymer_n_atoms, n_polymers, &
                                                  max_intersection_level, &
                                                  intersections, intersection_sets, intersection_levels, n_intersections)

      deallocate (polymer_atoms, polymer_n_atoms)
   end subroutine generate_polymer_intersections

   subroutine generate_intersections_from_atom_lists(atom_lists, n_atoms_list, n_sets, max_k_level, &
                                                   intersections, intersection_sets, intersection_levels, n_intersections)
      !! Generate k-way intersections from arbitrary atom lists (not tied to sys_geom)
      !! max_k_level limits the maximum intersection order to prevent combinatorial explosion
      use pic_logger, only: logger => global_logger
      use pic_io, only: to_char
      integer, intent(in) :: atom_lists(:, :)  !! (max_atoms, n_sets)
      integer, intent(in) :: n_atoms_list(:)   !! Number of atoms in each set
      integer, intent(in) :: n_sets            !! Number of sets (polymers)
      integer, intent(in) :: max_k_level       !! Maximum intersection level to compute
      integer, allocatable, intent(out) :: intersections(:, :)
      integer, allocatable, intent(out) :: intersection_sets(:, :)
      integer, allocatable, intent(out) :: intersection_levels(:)
      integer, intent(out) :: n_intersections

      integer :: max_intersections, max_atoms
      integer, allocatable :: temp_intersections(:, :), temp_sets(:, :), temp_levels(:)
      integer :: intersection_count, k, idx, actual_max_k
      integer, allocatable :: combination(:)

      if (n_sets < 2) then
         n_intersections = 0
         return
      end if

      ! Limit k-way intersections to min(max_k_level, n_sets)
      actual_max_k = min(max_k_level, n_sets)

      max_intersections = 2**n_sets - n_sets - 1
      max_atoms = maxval(n_atoms_list)

      allocate (temp_intersections(max_atoms, max_intersections))
      allocate (temp_sets(n_sets, max_intersections))
      allocate (temp_levels(max_intersections))
      temp_intersections = 0
      temp_sets = 0
      intersection_count = 0

      call logger%info("Generating k-way intersections (k=2 to "//to_char(actual_max_k)//")")

      ! Loop over intersection levels k from 2 to actual_max_k
      do k = 2, actual_max_k
         allocate (combination(k))
         call generate_k_way_intersections_from_lists(atom_lists, n_atoms_list, n_sets, k, &
                                                      combination, max_atoms, &
                                                      temp_intersections, temp_sets, temp_levels, intersection_count)
         deallocate (combination)
      end do

      n_intersections = intersection_count

      ! Allocate output arrays
      if (n_intersections > 0) then
         allocate (intersections(max_atoms, n_intersections))
         allocate (intersection_sets(n_sets, n_intersections))
         allocate (intersection_levels(n_intersections))
         intersections = temp_intersections(1:max_atoms, 1:n_intersections)
         intersection_sets = temp_sets(1:n_sets, 1:n_intersections)
         intersection_levels = temp_levels(1:n_intersections)

         call logger%info("Generated "//to_char(n_intersections)//" total intersections:")
         do k = 2, actual_max_k
            idx = count(intersection_levels == k)
            if (idx > 0) then
               call logger%info("  "//to_char(idx)//" intersections at level "//to_char(k))
            end if
         end do
      else
         allocate (intersections(1, 0))
         allocate (intersection_sets(1, 0))
         allocate (intersection_levels(0))
         call logger%info("No intersections found")
      end if

      deallocate (temp_intersections, temp_sets, temp_levels)
   end subroutine generate_intersections_from_atom_lists

   subroutine generate_k_way_intersections_from_lists(atom_lists, n_atoms_list, n_sets, k, combination, max_atoms, &
                                                      temp_intersections, temp_sets, temp_levels, intersection_count)
      !! Generate all k-way intersections from atom lists
      integer, intent(in) :: atom_lists(:, :), n_atoms_list(:), n_sets, k, max_atoms
      integer, intent(inout) :: combination(:)
      integer, intent(inout) :: temp_intersections(:, :), temp_sets(:, :), temp_levels(:)
      integer, intent(inout) :: intersection_count

      integer, allocatable :: current_intersection(:), temp_intersection(:)
      integer :: current_n_intersect, temp_n_intersect
      integer :: i, j, set_idx
      logical :: has_intersection

      call next_combination_init(combination, k)

      do
         ! Compute intersection of all sets in this combination
         has_intersection = .false.

         ! Start with the first set in the combination
         set_idx = combination(1)
         allocate (current_intersection(n_atoms_list(set_idx)))
         current_intersection = atom_lists(1:n_atoms_list(set_idx), set_idx)
         current_n_intersect = n_atoms_list(set_idx)

         ! Intersect with each subsequent set
         do i = 2, k
            set_idx = combination(i)
            allocate (temp_intersection(current_n_intersect))

            ! Find intersection
            temp_n_intersect = 0
            do j = 1, current_n_intersect
               if (any(atom_lists(1:n_atoms_list(set_idx), set_idx) == current_intersection(j))) then
                  temp_n_intersect = temp_n_intersect + 1
                  temp_intersection(temp_n_intersect) = current_intersection(j)
               end if
            end do

            ! Update current intersection
            deallocate (current_intersection)
            if (temp_n_intersect > 0) then
               allocate (current_intersection(temp_n_intersect))
               current_intersection = temp_intersection(1:temp_n_intersect)
               current_n_intersect = temp_n_intersect
               has_intersection = .true.
            else
               has_intersection = .false.
            end if
            deallocate (temp_intersection)

            if (.not. has_intersection) exit
         end do

         ! Store if we found an intersection
         if (has_intersection .and. current_n_intersect > 0) then
            intersection_count = intersection_count + 1
            temp_intersections(1:current_n_intersect, intersection_count) = current_intersection(1:current_n_intersect)
            temp_sets(1:k, intersection_count) = combination(1:k)
            temp_levels(intersection_count) = k
         end if

         if (allocated(current_intersection)) deallocate (current_intersection)

         ! Get next combination
         if (.not. next_combination(combination, k, n_sets)) exit
      end do
   end subroutine generate_k_way_intersections_from_lists

   subroutine gmbe_enumerate_pie_terms(sys_geom, primaries, n_primaries, polymer_level, max_k_level, &
                                       pie_atom_sets, pie_coefficients, n_pie_terms)
      !! Enumerate all unique intersections via DFS and accumulate PIE coefficients
      !! This implements the GMBE(N) algorithm with inclusion-exclusion principle
      !!
      !! Algorithm:
      !! 1. For each primary i, start DFS with clique=[i]
      !! 2. Recursively grow cliques by adding overlapping primaries
      !! 3. For each clique of size k, compute intersection and add PIE coefficient:
      !!    coefficient = (+1) if k odd, (-1) if k even
      !! 4. Accumulate coefficients for each unique atom set
      use mqc_physical_fragment, only: system_geometry_t

      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: primaries(:, :)   !! Primary polymers (n_primaries, polymer_level)
      integer, intent(in) :: n_primaries       !! Number of primary polymers
      integer, intent(in) :: polymer_level     !! Level of primaries (1=monomers, 2=dimers, etc.)
      integer, intent(in) :: max_k_level       !! Maximum clique size (intersection depth limit)
      integer, allocatable, intent(out) :: pie_atom_sets(:, :)  !! Unique atom sets (max_atoms, n_terms)
      integer, allocatable, intent(out) :: pie_coefficients(:)  !! PIE coefficient for each term
      integer, intent(out) :: n_pie_terms      !! Number of unique PIE terms

      ! Temporary storage for PIE terms (allocate generously)
      integer, parameter :: MAX_PIE_TERMS = 100000  ! Adjust if needed
      integer, allocatable :: temp_atom_sets(:, :)
      integer, allocatable :: temp_coefficients(:)
      integer, allocatable :: primary_atoms(:, :)    !! Precomputed atom lists for each primary
      integer, allocatable :: primary_n_atoms(:)     !! Atom counts for each primary
      integer, allocatable :: clique(:)              !! Current clique being built
      integer, allocatable :: current_atoms(:)       !! Current intersection atoms
      integer, allocatable :: candidates(:)          !! Candidate primaries to add
      integer :: i, j, max_atoms, n_candidates

      call logger%info("Enumerating GMBE PIE terms via DFS...")

      ! Find maximum atoms in any fragment
      max_atoms = sys_geom%total_atoms

      ! Allocate temporary storage
      allocate (temp_atom_sets(max_atoms, MAX_PIE_TERMS))
      allocate (temp_coefficients(MAX_PIE_TERMS))
      temp_atom_sets = -1
      temp_coefficients = 0
      n_pie_terms = 0

      ! Precompute atom lists for all primaries
      allocate (primary_atoms(max_atoms, n_primaries))
      allocate (primary_n_atoms(n_primaries))
      primary_atoms = -1

      do i = 1, n_primaries
         block
            integer, allocatable :: atom_list(:)
            integer :: n_atoms
            call compute_polymer_atoms(sys_geom, primaries(i, :), polymer_level, atom_list, n_atoms)
            primary_n_atoms(i) = n_atoms
            primary_atoms(1:n_atoms, i) = atom_list(1:n_atoms)
            deallocate (atom_list)
         end block
      end do

      ! Allocate work arrays
      allocate (clique(max_k_level))
      allocate (current_atoms(max_atoms))
      allocate (candidates(n_primaries))

      ! Start DFS from each primary
      do i = 1, n_primaries
         ! Initial clique: just primary i
         clique(1) = i
         current_atoms(1:primary_n_atoms(i)) = primary_atoms(1:primary_n_atoms(i), i)

         ! Candidates: all primaries after i (to avoid duplicates)
         n_candidates = n_primaries - i
         if (n_candidates > 0) then
            candidates(1:n_candidates) = [(i + j, j=1, n_candidates)]
         end if

         ! DFS from this primary
         call dfs_pie_accumulate(primary_atoms, primary_n_atoms, n_primaries, max_atoms, &
                                 clique, 1, current_atoms(1:primary_n_atoms(i)), primary_n_atoms(i), &
                                 candidates, n_candidates, max_k_level, &
                                 temp_atom_sets, temp_coefficients, n_pie_terms, MAX_PIE_TERMS)
      end do

      ! Copy to output arrays
      if (n_pie_terms > 0) then
         allocate (pie_atom_sets(max_atoms, n_pie_terms))
         allocate (pie_coefficients(n_pie_terms))
         pie_atom_sets = temp_atom_sets(:, 1:n_pie_terms)
         pie_coefficients = temp_coefficients(1:n_pie_terms)
      end if

      call logger%info("Generated "//to_char(n_pie_terms)//" unique PIE terms")

      ! Cleanup
      deallocate (temp_atom_sets, temp_coefficients, primary_atoms, primary_n_atoms)
      deallocate (clique, current_atoms, candidates)

   end subroutine gmbe_enumerate_pie_terms

   recursive subroutine dfs_pie_accumulate(primary_atoms, primary_n_atoms, n_primaries, max_atoms, &
                                           clique, clique_size, current_atoms, n_current_atoms, &
                                           candidates, n_candidates, max_k_level, &
                                           atom_sets, coefficients, n_terms, max_terms)
      !! DFS helper: accumulate PIE coefficients for intersections
      integer, intent(in) :: primary_atoms(:, :)    !! Precomputed atom lists
      integer, intent(in) :: primary_n_atoms(:)     !! Atom counts
      integer, intent(in) :: n_primaries, max_atoms
      integer, intent(in) :: clique(:)              !! Current clique
      integer, intent(in) :: clique_size            !! Size of current clique
      integer, intent(in) :: current_atoms(:)       !! Atoms in current intersection
      integer, intent(in) :: n_current_atoms        !! Number of atoms in intersection
      integer, intent(in) :: candidates(:)          !! Candidate primaries
      integer, intent(in) :: n_candidates
      integer, intent(in) :: max_k_level
      integer, intent(inout) :: atom_sets(:, :)
      integer, intent(inout) :: coefficients(:)
      integer, intent(inout) :: n_terms
      integer, intent(in) :: max_terms

      integer :: sign, term_idx, i, candidate_idx
      integer, allocatable :: new_atoms(:), new_candidates(:)
      integer :: n_new_atoms, n_new_candidates
      logical :: found

      ! Skip empty intersections
      if (n_current_atoms == 0) return

      ! Compute PIE sign: (+1) for odd clique size, (-1) for even
      sign = merge(1, -1, mod(clique_size, 2) == 1)

      ! Find or create entry for this atom set
      found = .false.
      do i = 1, n_terms
         if (atom_sets_equal(atom_sets(:, i), current_atoms, n_current_atoms)) then
            coefficients(i) = coefficients(i) + sign
            found = .true.
            exit
         end if
      end do

      if (.not. found) then
         ! New atom set
         if (n_terms >= max_terms) then
            call logger%error("Exceeded maximum PIE terms ("//to_char(max_terms)//")")
            return
         end if
         n_terms = n_terms + 1
         atom_sets(1:n_current_atoms, n_terms) = current_atoms(1:n_current_atoms)
         atom_sets(n_current_atoms + 1:, n_terms) = -1
         coefficients(n_terms) = sign
      end if

      ! Stop if we've reached maximum clique size
      if (clique_size >= max_k_level) return
      if (n_candidates == 0) return

      ! Try adding each candidate to the clique
      allocate (new_atoms(max_atoms))
      allocate (new_candidates(n_primaries))

      do i = 1, n_candidates
         candidate_idx = candidates(i)

         ! Compute intersection with this candidate
         call intersect_atom_lists(current_atoms, n_current_atoms, &
                                   primary_atoms(:, candidate_idx), primary_n_atoms(candidate_idx), &
                                   new_atoms, n_new_atoms)

         ! Skip if no intersection
         if (n_new_atoms == 0) cycle

         ! New candidates: must come after this one and overlap with new_atoms
         n_new_candidates = 0
         do term_idx = i + 1, n_candidates
            block
               integer :: test_candidate, test_n_intersect
               integer, allocatable :: test_intersect(:)
               test_candidate = candidates(term_idx)

               allocate (test_intersect(max_atoms))
               call intersect_atom_lists(new_atoms, n_new_atoms, &
                                         primary_atoms(:, test_candidate), primary_n_atoms(test_candidate), &
                                         test_intersect, test_n_intersect)

               if (test_n_intersect > 0) then
                  n_new_candidates = n_new_candidates + 1
                  new_candidates(n_new_candidates) = test_candidate
               end if
               deallocate (test_intersect)
            end block
         end do

         ! Recurse
         block
            integer :: new_clique(clique_size + 1)
            new_clique(1:clique_size) = clique(1:clique_size)
            new_clique(clique_size + 1) = candidate_idx

            call dfs_pie_accumulate(primary_atoms, primary_n_atoms, n_primaries, max_atoms, &
                                    new_clique, clique_size + 1, new_atoms, n_new_atoms, &
                                    new_candidates, n_new_candidates, max_k_level, &
                                    atom_sets, coefficients, n_terms, max_terms)
         end block
      end do

      deallocate (new_atoms, new_candidates)

   end subroutine dfs_pie_accumulate

   function atom_sets_equal(set1, set2, n_atoms) result(equal)
      !! Check if two atom sets are equal (assuming sorted)
      integer, intent(in) :: set1(:), set2(:)
      integer, intent(in) :: n_atoms
      logical :: equal
      integer :: i

      equal = .true.
      do i = 1, n_atoms
         if (set1(i) /= set2(i)) then
            equal = .false.
            return
         end if
      end do
   end function atom_sets_equal

   subroutine intersect_atom_lists(atoms1, n1, atoms2, n2, intersection, n_intersect)
      !! Compute intersection of two atom lists
      integer, intent(in) :: atoms1(:), n1, atoms2(:), n2
      integer, intent(out) :: intersection(:)
      integer, intent(out) :: n_intersect
      integer :: i, j

      n_intersect = 0
      do i = 1, n1
         if (atoms1(i) < 0) cycle
         do j = 1, n2
            if (atoms2(j) < 0) cycle
            if (atoms1(i) == atoms2(j)) then
               n_intersect = n_intersect + 1
               intersection(n_intersect) = atoms1(i)
               exit
            end if
         end do
      end do
   end subroutine intersect_atom_lists

end module mqc_frag_utils
