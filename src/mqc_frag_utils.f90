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

end module mqc_frag_utils
