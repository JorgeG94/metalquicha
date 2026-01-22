!! Hash-based lookup table for fast fragment index retrieval
module mqc_fragment_lookup
   !! Provides O(1) hash table for mapping monomer combinations to fragment indices
   use pic_types, only: int32, int64, dp
   use pic_sorting, only: sort
   use pic_hash_32bit, only: fnv_1a_hash
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: fragment_lookup_t  !! Hash-based lookup table type

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

   pure subroutine fragment_lookup_init(this, estimated_entries)
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

   subroutine fragment_lookup_insert(this, monomers, n, fragment_idx, error)
      !! Insert a monomer combination -> fragment index mapping
      class(fragment_lookup_t), intent(inout) :: this
      integer, intent(in) :: monomers(:), n
      integer(int64), intent(in) :: fragment_idx
      type(error_t), intent(out), optional :: error

      integer(int32) :: hash_val
      integer :: bucket
      type(hash_entry_t), pointer :: new_entry
      integer, allocatable :: sorted_key(:)

      if (.not. this%initialized) then
         if (present(error)) then
            call error%set(ERROR_VALIDATION, "Hash table not initialized")
         end if
         return
      end if

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

   pure subroutine fragment_lookup_destroy(this)
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

end module mqc_fragment_lookup
