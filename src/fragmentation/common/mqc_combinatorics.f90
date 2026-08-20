!! Combinatorial mathematics utilities for fragment generation
module mqc_combinatorics
   !! Provides pure combinatorial functions for generating molecular fragments
   !! including binomial coefficients, combinations, and fragment counting
   use pic_types, only: default_int, int32, int64
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_math_utils, only: binomial
   implicit none
   private

   public :: fragment_size_of      !! How many monomers a polymer row names
   public :: vmfc_subset_key       !! Counterpoise subset key: chosen real, rest ghosted
   public :: is_auxiliary_row      !! A ghosted row: subtracted, never summed
   public :: real_count_of         !! Real (non-ghosted) monomers in a row
   public :: binomial              !! Binomial coefficient calculation
   public :: get_nfrags            !! Calculate total number of fragments
   public :: create_monomer_list   !! Generate sequential monomer indices
   public :: generate_fragment_list  !! Generate all fragments up to max level
   public :: combine               !! Generate all combinations of size r
   public :: get_next_combination  !! Generate next combination in sequence
   public :: next_combination_init  !! Initialize combination to [1,2,...,k]
   public :: next_combination      !! Generate next combination (alternate interface)
   public :: print_combos          !! Debug utility to print combinations
   public :: calculate_fragment_distances  !! Calculate minimal distances for all fragments

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

   pure subroutine create_monomer_list(monomers)
      !! Generate a list of monomer indices from 1 to N
      integer(default_int), allocatable, intent(inout) :: monomers(:)
      integer(default_int) :: i, length

      length = size(monomers, 1)

      do i = 1, length
         monomers(i) = i
      end do

   end subroutine create_monomer_list

   pure function is_auxiliary_row(row) result(aux)
      !! Whether a row exists only to be subtracted, not to be summed
      !!
      !! A counterpoise expansion computes monomer A in the basis of the pair
      !! AB. That energy belongs inside the pair's correction and nowhere else:
      !!
      !!     E = sum_i E_i(i)  +  sum_ij [ E_ij - E_i(ij) - E_j(ij) ]
      !!
      !! The one-body term uses each monomer in its *own* basis, so the ghosted
      !! rows are auxiliary -- adding their deltas to the total as well would
      !! count them twice. A negative entry is what marks them, which is the
      !! same sign the fragment builder and the subset key already read.
      integer(default_int), intent(in) :: row(:)
      logical :: aux

      aux = any(row < 0)
   end function is_auxiliary_row

   pure function real_count_of(row) result(n)
      !! How many of a row's monomers are real rather than ghosted
      !!
      !! The subset recursion works over these: [1,-2] contains one real
      !! monomer, so it has no proper subsets and its delta is its energy. Using
      !! the full row size instead would send the recursion looking for subsets
      !! that were never generated.
      integer(default_int), intent(in) :: row(:)
      integer(default_int) :: n

      n = count(row > 0)
   end function real_count_of

   pure subroutine vmfc_subset_key(fragment, n, chosen, k, key)
      !! The subset key a counterpoise-corrected expansion looks up
      !!
      !! Ordinary MBE subtracts the subset {A} from the pair {A,B}. VMFC
      !! subtracts {A in the basis of AB} instead -- the same monomers, solved
      !! in the parent's basis -- so the superposition error that inflates the
      !! parent stands on both sides of the difference and cancels rather than
      !! surviving into the total.
      !!
      !! So the key is the chosen monomers positive and *everything else in the
      !! parent* negative. For the pair [1,2] choosing [1], that is [1,-2]. The
      !! recursion, the hash table and the dispatch are untouched: a signed key
      !! is just another key, and it sorts and hashes like any other.
      integer, intent(in) :: fragment(:)   !! The parent's monomers, all positive
      integer, intent(in) :: n             !! How many of them
      integer, intent(in) :: chosen(:)     !! Positions within `fragment`, size k
      integer, intent(in) :: k
      integer, intent(out) :: key(:)       !! Size n: k real, then n-k ghosted

      integer :: i, j, next
      logical :: taken(n)

      taken = .false.
      do i = 1, k
         taken(chosen(i)) = .true.
         key(i) = fragment(chosen(i))
      end do

      next = k
      do j = 1, n
         if (.not. taken(j)) then
            next = next + 1
            key(next) = -fragment(j)
         end if
      end do
   end subroutine vmfc_subset_key

   pure function fragment_size_of(row) result(n)
      !! How many monomers a polymer row names, padding excluded
      !!
      !! Rows are zero-padded to the widest fragment, so the size is the count
      !! of non-zero entries. Written `/= 0` rather than `> 0` because a
      !! negative entry is a monomer present as *ghost centres* -- its atoms and
      !! basis functions without its nucleus -- which is a monomer the fragment
      !! contains and must be sized for. `> 0` would drop it silently and hand
      !! the SCF a smaller molecule than the row asked for.
      !!
      !! Here rather than repeated at each site because there were seventeen of
      !! them, and one missed would be a wrong energy rather than a failure.
      integer(default_int), intent(in) :: row(:)
      integer(default_int) :: n

      n = count(row /= 0)
   end function fragment_size_of

   recursive subroutine generate_fragment_list(monomers, max_level, polymers, count)
      !! Generate all possible fragments (combinations of monomers) up to max_level
      !! Uses int64 for count to handle large numbers of fragments that overflow int32.
      integer(default_int), intent(in) ::  max_level
      integer(default_int), intent(in) :: monomers(:)
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
      integer(default_int), intent(in) ::  n, r, index, i
      integer(default_int), intent(in) :: arr(:)
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
      integer(default_int), intent(in) ::  max_len
      integer(default_int), intent(in) :: out_array(:, :)
      integer(int64), intent(in) :: count
      integer(int64) :: i
      integer(default_int) :: j

      character(len=MAX_LINE_LENGTH) :: line

      ! Assembled in `line` and emitted once per combination. The non-advancing
      ! writes this replaces built the row directly on the unit, which a logger
      ! cannot do: a log record is a whole line or it is nothing.
      do i = 1_int64, count
         line = ""
         do j = 1, max_len
            if (out_array(i, j) == 0) exit
            line = trim(line)//to_char(out_array(i, j))
            if (j < max_len .and. out_array(i, j + 1) /= 0) then
               line = trim(line)//":"
            end if
         end do
         call logger%info(trim(line))
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

   subroutine calculate_fragment_distances(polymers, fragment_count, sys_geom, distances)
      !! Calculate minimal atomic distance for each fragment
      !! For monomers (1-body), distance is 0.0
      !! For n-mers (n >= 2), distance is the minimum distance between atoms
      !! in different constituent monomers
      use pic_types, only: dp
      use mqc_physical_fragment, only: system_geometry_t, to_angstrom
      integer(default_int), intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: fragment_count
      type(system_geometry_t), intent(in) :: sys_geom
      real(dp), intent(out) :: distances(:)

      integer(int64) :: ifrag
      integer :: fragment_size, i, j, iatom, jatom
      integer :: mon_i, mon_j
      integer :: atom_start_i, atom_end_i, atom_start_j, atom_end_j
      integer :: k
      real(dp) :: dist, min_dist
      real(dp) :: dx, dy, dz
      logical :: is_variable_size

      ! Check if we have variable-sized fragments
      is_variable_size = allocated(sys_geom%fragment_sizes)

      do ifrag = 1_int64, fragment_count
         fragment_size = fragment_size_of(polymers(ifrag, :))

         if (fragment_size == 1) then
            ! Monomers have distance 0
            distances(ifrag) = 0.0_dp
         else
            ! For n-mers, calculate minimal distance between atoms in different monomers
            min_dist = huge(1.0_dp)

            ! Loop over all pairs of monomers in this fragment
            do i = 1, fragment_size - 1
               mon_i = polymers(ifrag, i)
               do j = i + 1, fragment_size
                  mon_j = polymers(ifrag, j)

                  if (is_variable_size) then
                     ! Variable-sized fragments: use fragment_atoms to get atom indices
                     ! Count atoms in this fragment
                     do iatom = 1, sys_geom%fragment_sizes(mon_i)
                        atom_start_i = sys_geom%fragment_atoms(iatom, mon_i) + 1  ! Convert to 1-indexed
                        do jatom = 1, sys_geom%fragment_sizes(mon_j)
                           atom_start_j = sys_geom%fragment_atoms(jatom, mon_j) + 1  ! Convert to 1-indexed

                           ! Calculate distance
                           dx = sys_geom%coordinates(1, atom_start_i) - sys_geom%coordinates(1, atom_start_j)
                           dy = sys_geom%coordinates(2, atom_start_i) - sys_geom%coordinates(2, atom_start_j)
                           dz = sys_geom%coordinates(3, atom_start_i) - sys_geom%coordinates(3, atom_start_j)
                           dist = sqrt(dx*dx + dy*dy + dz*dz)

                           if (dist < min_dist) min_dist = dist
                        end do
                     end do
                  else
                     ! Fixed-sized monomers: calculate atom range directly
                     atom_start_i = (mon_i - 1)*sys_geom%atoms_per_monomer + 1
                     atom_end_i = mon_i*sys_geom%atoms_per_monomer

                     atom_start_j = (mon_j - 1)*sys_geom%atoms_per_monomer + 1
                     atom_end_j = mon_j*sys_geom%atoms_per_monomer

                     ! Loop over all atoms in monomer i
                     do iatom = atom_start_i, atom_end_i
                        ! Loop over all atoms in monomer j
                        do jatom = atom_start_j, atom_end_j
                           ! Calculate distance (coordinates are in Bohr)
                           dx = sys_geom%coordinates(1, iatom) - sys_geom%coordinates(1, jatom)
                           dy = sys_geom%coordinates(2, iatom) - sys_geom%coordinates(2, jatom)
                           dz = sys_geom%coordinates(3, iatom) - sys_geom%coordinates(3, jatom)
                           dist = sqrt(dx*dx + dy*dy + dz*dz)

                           if (dist < min_dist) min_dist = dist
                        end do
                     end do
                  end if
               end do
            end do

            ! Convert from Bohr to Angstrom
            distances(ifrag) = to_angstrom(min_dist)
         end if
      end do

   end subroutine calculate_fragment_distances

end module mqc_combinatorics
