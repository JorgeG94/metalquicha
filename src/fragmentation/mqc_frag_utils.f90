!! Fragment generation and manipulation utilities
module mqc_frag_utils
   !! Provides combinatorial functions and algorithms for generating molecular
   !! fragments, managing fragment lists, and performing many-body expansion calculations.
   use pic_types, only: default_int, dp
   use pic_logger, only: logger => global_logger
   implicit none
   private
   public :: binomial              !! Binomial coefficient calculation
   public :: create_monomer_list   !! Generate sequential monomer indices
   public :: generate_fragment_list  !! Generate all fragments up to max level
   public :: get_nfrags            !! Calculate total number of fragments
   public :: next_combination      !! Generate next combination in sequence
   public :: find_fragment_index   !! Locate fragment by composition

contains

   pure function get_nfrags(n_monomers, max_level) result(n_expected_fragments)
      !! Calculate total number of fragments for given system size and max level
      !!
      !! Computes the sum of binomial coefficients C(n,k) for k=1 to max_level,
      !! representing all possible fragments from monomers to max_level-mers.
      integer(default_int), intent(in) :: n_monomers  !! Number of monomers in system
      integer(default_int), intent(in) :: max_level   !! Maximum fragment size
      integer(default_int) :: n_expected_fragments     !! Total fragment count
      integer(default_int) :: i  !! Loop counter

      n_expected_fragments = 0
      do i = 1, max_level
         n_expected_fragments = n_expected_fragments + binomial(n_monomers, i)
      end do
   end function get_nfrags

   pure function binomial(n, r) result(c)
      !! Compute binomial coefficient C(n,r) = n! / (r! * (n-r)!)
      !!
      !! Calculates "n choose r" using iterative algorithm to avoid
      !! factorial overflow for large numbers.
      integer(default_int), intent(in) :: n  !! Total number of items
      integer(default_int), intent(in) :: r  !! Number of items to choose
      integer(default_int) :: c              !! Binomial coefficient result
      integer(default_int) :: i              !! Loop counter

      if (r == 0 .or. r == n) then
         c = 1
      else if (r > n) then
         c = 0
      else
         c = 1
         do i = 1, r
            c = c*(n - i + 1)/i
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
      integer(default_int), intent(in) :: monomers(:), max_level
      integer(default_int), intent(inout) :: polymers(:, :)
      integer(default_int), intent(inout) :: count
      integer(default_int) :: r, n

      n = size(monomers, 1)
      do r = 2, max_level
         call combine(monomers, n, r, polymers, count)
      end do
   end subroutine generate_fragment_list

   recursive subroutine combine(arr, n, r, out_array, count)
      !! Generate all combinations of size r from array arr of size n
      integer(default_int), intent(in) :: arr(:)
      integer(default_int), intent(in) :: n, r
      integer(default_int), intent(inout) :: out_array(:, :)
      integer(default_int), intent(inout) :: count
      integer(default_int) :: data(r)
      call combine_util(arr, n, r, 1, data, 1, out_array, count)
   end subroutine combine

   recursive subroutine combine_util(arr, n, r, index, data, i, out_array, count)
      !! Utility for generating combinations recursively
      integer(default_int), intent(in) :: arr(:), n, r, index, i
      integer(default_int), intent(inout) :: data(:), out_array(:, :)
      integer(default_int), intent(inout) :: count
      integer(default_int) :: j

      if (index > r) then
         count = count + 1
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
      integer(default_int), intent(in) :: out_array(:, :), count, max_len
      integer(default_int) :: i, j

      do i = 1, count
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

   function next_combination(indices, k, n) result(has_next)
      !! Generate next combination (updates indices in place)
      !! Returns .true. if there's a next combination, .false. if we're done
      integer, intent(inout) :: indices(:)
      integer, intent(in) :: k, n
      logical :: has_next
      integer :: i

      has_next = .true.

      ! Find rightmost index that can be incremented
      i = k
      do while (i >= 1)
         if (indices(i) < n - k + i) then
            indices(i) = indices(i) + 1
            ! Reset all indices to the right
            do while (i < k)
               i = i + 1
               indices(i) = indices(i - 1) + 1
            end do
            return
         end if
         i = i - 1
      end do

      ! No more combinations
      has_next = .false.

   end function next_combination

   function find_fragment_index(target_monomers, polymers, fragment_count, expected_size) result(idx)
      !! Find the fragment index that contains exactly the target monomers
      integer, intent(in) :: target_monomers(:), polymers(:, :), fragment_count, expected_size
      integer :: idx

      integer :: i, j, fragment_size
      logical :: match

      idx = -1

      do i = 1, fragment_count
         fragment_size = count(polymers(i, :) > 0)

         ! Check if this fragment has the right size
         if (fragment_size /= expected_size) cycle

         ! Check if all target monomers are in this fragment
         match = .true.
         do j = 1, expected_size
            if (.not. any(polymers(i, 1:fragment_size) == target_monomers(j))) then
               match = .false.
               exit
            end if
         end do

         if (match) then
            idx = i
            return
         end if
      end do

      ! If we get here, we didn't find the fragment
      block
         character(len=256) :: monomers_str
         integer :: k
         write (monomers_str, '(*(i0,1x))') (target_monomers(k), k=1, size(target_monomers))
         call logger%error("Could not find fragment with monomers: "//trim(monomers_str))
      end block
      error stop "Fragment not found in find_fragment_index"

   end function find_fragment_index

end module mqc_frag_utils
