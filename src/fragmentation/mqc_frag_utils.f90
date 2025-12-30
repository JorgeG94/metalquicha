!! Fragment generation and manipulation utilities
module mqc_frag_utils
   !! Provides combinatorial functions and algorithms for generating molecular
   !! fragments, managing fragment lists, and performing many-body expansion calculations.
   use pic_types, only: default_int, int32, int64, dp
   use pic_logger, only: logger => global_logger
   use pic_sorting, only: sort
   use pic_hash_32bit, only: fnv_1a_hash
   implicit none
   private
   public :: binomial              !! Binomial coefficient calculation
   public :: create_monomer_list   !! Generate sequential monomer indices
   public :: generate_fragment_list  !! Generate all fragments up to max level
   public :: get_nfrags            !! Calculate total number of fragments
   public :: get_next_combination      !! Generate next combination in sequence
   public :: fragment_lookup_t     !! Hash-based lookup table for fast fragment index retrieval
   public :: find_fragment_intersection  !! Find shared atoms between two fragments (for GMBE)
   public :: generate_intersections     !! Generate all pairwise intersections for GMBE

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

   subroutine generate_intersections(sys_geom, monomers, polymers, n_monomers, &
                                     intersections, intersection_pairs, n_intersections)
      !! Generate all pairwise intersections between fragments
      !!\n!! For a system with overlapping fragments, this identifies all pairs of fragments
      !! that share atoms and stores the intersection atom lists.
      !!\n!! Algorithm:
      !! - Loop through all pairs (i,j) where i < j
      !! - Call find_fragment_intersection for each pair
      !! - Store non-empty intersections and track which monomer pair created them
      use mqc_physical_fragment, only: system_geometry_t
      type(system_geometry_t), intent(in) :: sys_geom
      integer, intent(in) :: monomers(:)           !! Monomer indices
      integer, intent(inout) :: polymers(:, :)     !! Output: monomers stored here
      integer, intent(in) :: n_monomers            !! Number of monomers
      integer, allocatable, intent(out) :: intersections(:, :)  !! Intersection atom lists
      integer, allocatable, intent(out) :: intersection_pairs(:, :)  !! Which pair (i,j) created each intersection
      integer, intent(out) :: n_intersections      !! Number of intersections found

      ! Temporaries for storing intersections
      integer, allocatable :: temp_intersections(:, :)
      integer, allocatable :: temp_pairs(:, :)
      integer, allocatable :: temp_intersection(:)
      integer :: temp_n_intersect
      logical :: has_intersection
      integer :: i, j, pair_count, max_atoms, max_intersections
      integer :: mono_i, mono_j, size_i, size_j

      ! Store monomers in polymers array
      polymers(1:n_monomers, 1) = monomers(1:n_monomers)

      ! Count maximum possible intersections: C(n,2) = n*(n-1)/2
      max_intersections = (n_monomers*(n_monomers - 1))/2

      if (max_intersections == 0) then
         n_intersections = 0
         return
      end if

      ! Find maximum atoms in any fragment for allocation
      max_atoms = maxval(sys_geom%fragment_sizes(1:n_monomers))

      ! Allocate temporary arrays
      allocate (temp_intersections(max_atoms, max_intersections))
      allocate (temp_pairs(2, max_intersections))
      temp_intersections = 0
      pair_count = 0

      ! Loop through all pairs (i, j) where i < j
      do i = 1, n_monomers - 1
         do j = i + 1, n_monomers
            ! Get fragment atoms for this pair
            mono_i = monomers(i)
            mono_j = monomers(j)
            size_i = sys_geom%fragment_sizes(mono_i)
            size_j = sys_geom%fragment_sizes(mono_j)

            ! Find intersection
            has_intersection = find_fragment_intersection( &
                               sys_geom%fragment_atoms(1:size_i, mono_i), size_i, &
                               sys_geom%fragment_atoms(1:size_j, mono_j), size_j, &
                               temp_intersection, temp_n_intersect)

            if (has_intersection) then
               pair_count = pair_count + 1
               temp_pairs(1, pair_count) = mono_i
               temp_pairs(2, pair_count) = mono_j
               temp_intersections(1:temp_n_intersect, pair_count) = temp_intersection

               ! Warn if mixing charged monomers (intersection assumed neutral)
               if (allocated(sys_geom%fragment_charges)) then
                  if (sys_geom%fragment_charges(mono_i) /= 0 .and. &
                      sys_geom%fragment_charges(mono_j) /= 0) then
                     block
                        character(len=20) :: str_i, str_j
                        write (str_i, '(I0)') mono_i
                        write (str_j, '(I0)') mono_j
                        call logger%warning("generate_intersections", &
                                            "Creating intersection between two charged monomers "// &
                                            "(fragments "//trim(str_i)//" and "//trim(str_j)//"). "// &
                                            "Intersection fragment will be NEUTRAL (charge=0), which may be unphysical.")
                     end block
                  end if
               end if

               deallocate (temp_intersection)
            end if
         end do
      end do

      n_intersections = pair_count

      ! Allocate output arrays
      if (n_intersections > 0) then
         allocate (intersections(max_atoms, n_intersections))
         allocate (intersection_pairs(2, n_intersections))
         intersections = temp_intersections(1:max_atoms, 1:n_intersections)
         intersection_pairs = temp_pairs(1:2, 1:n_intersections)
      end if

      deallocate (temp_intersections, temp_pairs)

   end subroutine generate_intersections

end module mqc_frag_utils
