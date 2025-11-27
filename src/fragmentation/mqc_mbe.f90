module mqc_frag_utils
   use pic_types, only: default_int, dp
   implicit none
   private
   public :: binomial
   public :: create_monomer_list
   public :: generate_fragment_list
   public :: get_nfrags

contains

   pure function get_nfrags(n_monomers, max_level) result(n_expected_fragments)
      integer(default_int), intent(in) :: n_monomers, max_level
      integer(default_int) :: n_expected_fragments
      integer(default_int) :: i

      n_expected_fragments = 0
      do i = 1, max_level
         n_expected_fragments = n_expected_fragments + binomial(n_monomers, i)
      end do
   end function get_nfrags

   pure function binomial(n, r) result(c)
      integer(default_int), intent(in) :: n, r
      integer(default_int) :: c
      integer(default_int) :: i

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

   subroutine create_monomer_list(monomers)
      integer(default_int), allocatable, intent(inout) :: monomers(:)
      integer(default_int) :: i, length

      length = size(monomers, 1)

      do i = 1, length
         monomers(i) = i
      end do

   end subroutine create_monomer_list

   recursive subroutine generate_fragment_list(monomers, max_level, polymers, count)
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
      integer(default_int), intent(in) :: arr(:)
      integer(default_int), intent(in) :: n, r
      integer(default_int), intent(inout) :: out_array(:, :)
      integer(default_int), intent(inout) :: count
      integer(default_int) :: data(r)
      call combine_util(arr, n, r, 1, data, 1, out_array, count)
   end subroutine combine

   recursive subroutine combine_util(arr, n, r, index, data, i, out_array, count)
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

end module mqc_frag_utils
