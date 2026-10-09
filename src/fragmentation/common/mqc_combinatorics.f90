!! Combinatorial mathematics utilities for fragment generation
module mqc_combinatorics
   !! Provides pure combinatorial functions for generating molecular fragments
   !! including binomial coefficients, combinations, and fragment counting
   use pic_types, only: default_int, int32, int64
   use pic_sorting, only: sort
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   use mqc_program_limits, only: MAX_LINE_LENGTH
   use mqc_math_utils, only: binomial
   implicit none
   private

   public :: fragment_size_of      !! How many monomers a polymer row names
   public :: vmfc_subset_key       !! Counterpoise subset key: chosen real, rest ghosted
   public :: vmfc_row_subset_key   !! The same, for a row that may already carry ghosts
   public :: is_auxiliary_row      !! A row subtracted by its parent, never summed
   public :: rows_match_counterpoise  !! Whether a term list is what a scheme produces
   public :: counterpoise_scheme_of  !! The scheme constant a deck's spelling names
   public :: counterpoise_scheme_name  !! The spelling of a scheme constant, for messages
   public :: COUNTERPOISE_NONE     !! Every subfragment in its own basis
   public :: COUNTERPOISE_VMFC     !! Every subfragment in its parent's basis
   public :: COUNTERPOISE_SSFC     !! Every subfragment in the whole cluster's basis
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

   integer, parameter :: COUNTERPOISE_NONE = 0
   integer, parameter :: COUNTERPOISE_VMFC = 1
   integer, parameter :: COUNTERPOISE_SSFC = 2

contains

   pure function get_nfrags(n_monomers, max_level) result(n_expected_fragments)
      !! Total number of fragments for a system size and a maximum level
      !!
      !! The sum of `C(n, k)` for `k = 1` to `max_level`.
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

   pure function counterpoise_scheme_of(name) result(scheme)
      !! The `COUNTERPOISE_*` constant for a `counterpoise` spelling
      !!
      !! `"vmfc"` and `"ssfc"` name their schemes; every other spelling,
      !! `"none"` included, gives `COUNTERPOISE_NONE`.
      ! `check_counterpoise_support` refuses every spelling it does not know
      ! before a term list is built, so the fall-through is not reached by a
      ! deck that got that far. It does not accept "ssfc" yet.
      character(len=*), intent(in) :: name
      integer :: scheme

      select case (trim(name))
      case ("vmfc")
         scheme = COUNTERPOISE_VMFC
      case ("ssfc")
         scheme = COUNTERPOISE_SSFC
      case default
         scheme = COUNTERPOISE_NONE
      end select
   end function counterpoise_scheme_of

   pure function counterpoise_scheme_name(scheme) result(name)
      !! The deck's spelling of a `COUNTERPOISE_*` constant
      integer, intent(in) :: scheme
      character(len=:), allocatable :: name

      select case (scheme)
      case (COUNTERPOISE_NONE)
         name = "none"
      case (COUNTERPOISE_VMFC)
         name = "vmfc"
      case (COUNTERPOISE_SSFC)
         name = "ssfc"
      case default
         name = "unknown"
      end select
   end function counterpoise_scheme_name

   pure function is_auxiliary_row(row, scheme) result(aux)
      !! Whether a row exists only to be subtracted, not to be summed
      !!
      !! A counterpoise expansion computes monomer A in a basis larger than its
      !! own. That energy belongs inside the correction of the term that
      !! subtracts it, and nowhere else; adding its delta to the total as well
      !! would count it twice. A negative entry marks a ghosted monomer, and
      !! what a ghost means for the sum depends on the scheme:
      !!
      !!     row                     VMFC        SSFC
      !!     all real                summed      summed
      !!     ghosted, 1 real         auxiliary   auxiliary
      !!     ghosted, 2+ real        auxiliary   summed
      !!
      !! Under VMFC every ghosted row is a subset solved in its parent's basis,
      !! so none is a term of the expansion. Under SSFC every term of order 2
      !! and above is itself a ghosted row, in the basis of the whole cluster;
      !! only the one-real-monomer rows are subtracted and not summed, since
      !! the one-body term is each monomer in its *own* basis. `COUNTERPOISE_NONE`
      !! has no ghosted rows to classify.
      integer(default_int), intent(in) :: row(:)
      integer, intent(in) :: scheme
         !! One of the `COUNTERPOISE_*` constants
      logical :: aux

      select case (scheme)
      case (COUNTERPOISE_VMFC)
         aux = any(row < 0)
      case (COUNTERPOISE_SSFC)
         aux = any(row < 0) .and. count(row > 0) == 1
      case default
         aux = .false.
      end select
   end function is_auxiliary_row

   pure function real_count_of(row) result(n)
      !! How many of a row's monomers are real rather than ghosted
      !!
      !! The subset recursion works over these: `[1,-2]` contains one real
      !! monomer, so it has no proper subsets and its delta is its energy.
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
      !! parent* negative. For the pair `[1,2]` choosing `[1]`, that is
      !! `[1,-2]`.
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

   pure subroutine vmfc_row_subset_key(row, chosen, k, key, key_len)
      !! The counterpoise subset key of a term-list row, ghosts included
      !!
      !! `chosen` picks among the row's *real* monomers, in the order they
      !! appear. The key keeps those real, ghosts the row's other real
      !! monomers, and keeps every ghost the row already carries. So the row
      !! `[1,2,-3]` choosing its first real monomer gives `[1,-2,-3]`: monomer 1
      !! in the basis of all three, which is the subset Valiron-Mayer subtracts.
      !! Member order is not significant; the lookup sorts keys.
      integer, intent(in) :: row(:)
         !! Zero-padded; positive entries real, negative ghosted
      integer, intent(in) :: chosen(:)   !! Positions among the real entries, size k
      integer, intent(in) :: k
      integer, intent(out) :: key(:)     !! At least as long as the row's non-zero entries
      integer, intent(out) :: key_len

      integer :: reals(size(row)), ghosts(size(row))
      integer :: n_real, n_ghost, i

      n_real = 0
      n_ghost = 0
      do i = 1, size(row)
         if (row(i) > 0) then
            n_real = n_real + 1
            reals(n_real) = row(i)
         else if (row(i) < 0) then
            n_ghost = n_ghost + 1
            ghosts(n_ghost) = row(i)
         end if
      end do

      call vmfc_subset_key(reals(1:n_real), n_real, chosen, k, key(1:n_real))
      key(n_real + 1:n_real + n_ghost) = ghosts(1:n_ghost)
      key_len = n_real + n_ghost
   end subroutine vmfc_row_subset_key

   pure function rows_match_counterpoise(polymers, n_rows, n_monomers, scheme) result(ok)
      !! Whether a term list is exactly what a counterpoise scheme produces
      !!
      !! The scheme is the caller's to name; the rows are checked against it by
      !! regeneration, rule by rule:
      !!
      !! `COUNTERPOISE_NONE`: no entry is negative.
      !!
      !! `COUNTERPOISE_VMFC`: the all-real rows are closed under taking one
      !! monomer away, every all-real row of size n >= 2 has all of its
      !! `2^n - 2` ghosted subsets as rows, and there are no other ghosted rows.
      !!
      !! `COUNTERPOISE_SSFC`: every ghosted row names all `n_monomers` monomers
      !! once, and no all-real row has 2 to `n_monomers - 1` monomers. Taking
      !! one real monomer from a ghosted or whole-system row and ghosting it
      !! gives another row of the list, down to the one-real-monomer rows,
      !! whose own-basis `[i]` row is present. Each own-basis monomer row has
      !! its full-basis row. A list with no ghosted row and no n-mer at all is
      !! level 1, which adds no rows under either scheme, and is accepted.
      !!
      !! A list with two equal rows is refused under every scheme but NONE.
      !! N = 2 is accepted by both VMFC and SSFC: the lists coincide.
      ! The scheme cannot be read off the rows. At L = N the SSFC top term is
      ! the whole system, whose complement is empty, so it is an all-real row
      ! exactly like VMFC's, and its ghosted rows are the ghosted subsets of
      ! that row. The VMFC closure rule is what tells the two lists apart: the
      ! SSFC list has none of the all-real n-mers a VMFC list holds beneath
      ! its top term.
      integer(default_int), intent(in) :: polymers(:, :)
         !! (rows, width), zero-padded; a negative entry is a ghosted monomer
      integer(int64), intent(in) :: n_rows  !! Rows of `polymers` in use
      integer, intent(in) :: n_monomers     !! Monomers in the whole system
      integer, intent(in) :: scheme         !! One of the `COUNTERPOISE_*` constants
      logical :: ok

      select case (scheme)
      case (COUNTERPOISE_NONE)
         ok = .not. any(polymers(1:n_rows, :) < 0)
      case (COUNTERPOISE_VMFC)
         ok = rows_match_vmfc(polymers(1:n_rows, :))
      case (COUNTERPOISE_SSFC)
         ok = rows_match_ssfc(polymers(1:n_rows, :), n_monomers)
      case default
         ok = .false.
      end select
   end function rows_match_counterpoise

   pure function rows_match_vmfc(rows) result(ok)
      !! The VMFC rule of `rows_match_counterpoise`
      !!
      !! Each ghosted row names a parent, its monomers taken as all real, and
      !! a real subset of it. Rows are distinct and a parent of n monomers has
      !! `2^n - 2` such subsets, so the ghosted rows are exactly the subsets
      !! when every one has a parent in the list and their number is the sum.
      integer(default_int), intent(in) :: rows(:, :)
      logical :: ok

      integer, allocatable :: keys(:, :)
      integer(int64), allocatable :: order(:)
      integer(int64) :: i, n_ghosted, n_expected
      integer :: n_real, p

      ok = .false.
      if (size(rows, 1) == 0) then
         ok = .true.
         return
      end if

      call canonical_rows(rows, keys)
      call order_keys(keys, order)
      if (has_equal_rows(keys, order)) return

      n_ghosted = 0_int64
      n_expected = 0_int64
      do i = 1_int64, size(rows, 1, kind=int64)
         n_real = count(rows(i, :) > 0)
         if (any(rows(i, :) < 0)) then
            if (n_real == 0) return
            if (.not. has_key(keys, order, absolute_key(keys(i, :)))) return
            n_ghosted = n_ghosted + 1_int64
         else if (n_real >= 2) then
            if (any(keys(i, 2:n_real) == keys(i, 1:n_real - 1))) return
            do p = 1, n_real
               if (.not. has_key(keys, order, without_entry(keys(i, :), p))) return
            end do
            n_expected = n_expected + 2_int64**n_real - 2_int64
         end if
      end do

      ok = n_ghosted == n_expected
   end function rows_match_vmfc

   pure function rows_match_ssfc(rows, n_monomers) result(ok)
      !! The SSFC rule of `rows_match_counterpoise`
      integer(default_int), intent(in) :: rows(:, :)
      integer, intent(in) :: n_monomers
      logical :: ok

      integer, allocatable :: keys(:, :), expected(:)
      integer(int64), allocatable :: order(:)
      integer(int64) :: i
      integer :: n_real, n_ghost, p
      logical :: has_ghosted

      ok = .false.
      if (size(rows, 1) == 0) then
         ok = .true.
         return
      end if

      call canonical_rows(rows, keys)
      call order_keys(keys, order)
      if (has_equal_rows(keys, order)) return
      has_ghosted = any(rows < 0)

      do i = 1_int64, size(rows, 1, kind=int64)
         n_real = count(rows(i, :) > 0)
         n_ghost = count(rows(i, :) < 0)
         if (n_ghost == 0 .and. n_real < 2) cycle      ! an own-basis monomer row

         ! The rest is a term in the whole cluster's basis, ghosted or, for the
         ! whole system at L = N, with nothing left to ghost.
         if (n_real == 0 .or. n_real + n_ghost /= n_monomers) return
         if (.not. names_whole_system(keys(i, :), n_monomers)) return

         if (n_real == 1) then
            ! Its own-basis row, which is where the correction starts from.
            expected = [keys(i, n_monomers), (0, p=1, size(keys, 2) - 1)]
            if (.not. has_key(keys, order, expected)) return
         else
            ! Every subset one monomer smaller is in the list.
            do p = 1, n_monomers
               if (keys(i, p) < 0) cycle
               if (.not. has_key(keys, order, ghosted_entry(keys(i, :), p))) return
            end do
         end if
      end do

      if (has_ghosted) then
         ! Every kept monomer has its full-basis row as well as its own.
         if (size(rows, 2) < n_monomers) return
         do i = 1_int64, size(rows, 1, kind=int64)
            if (count(rows(i, :) /= 0) /= 1 .or. any(rows(i, :) < 0)) cycle
            if (keys(i, 1) > n_monomers) return
            expected = full_basis_key(keys(i, 1), n_monomers, size(keys, 2))
            if (.not. has_key(keys, order, expected)) return
         end do
      end if

      ok = .true.
   end function rows_match_ssfc

   pure subroutine canonical_rows(rows, keys)
      !! Each row's non-zero entries in ascending order, zero-padded again
      !!
      !! Two rows name the same term exactly when their keys are equal. Ghosts
      !! sort before the real monomers, being negative.
      integer(default_int), intent(in) :: rows(:, :)
      integer, allocatable, intent(out) :: keys(:, :)

      integer, allocatable :: entries(:)
      integer(int64) :: i

      allocate (keys(size(rows, 1), size(rows, 2)))
      keys = 0
      do i = 1_int64, size(rows, 1, kind=int64)
         entries = pack(rows(i, :), rows(i, :) /= 0)
         call sort(entries)
         keys(i, 1:size(entries)) = entries
      end do
   end subroutine canonical_rows

   pure function key_before(a, b) result(before)
      !! Whether key `a` sorts strictly before key `b`, entry by entry
      integer, intent(in) :: a(:), b(:)
      logical :: before

      integer :: j

      before = .false.
      do j = 1, size(a)
         if (a(j) /= b(j)) then
            before = a(j) < b(j)
            return
         end if
      end do
   end function key_before

   pure subroutine order_keys(keys, order)
      !! The permutation that sorts the rows of `keys`, by merging runs
      integer, intent(in) :: keys(:, :)
      integer(int64), allocatable, intent(out) :: order(:)

      integer(int64), allocatable :: merged(:)
      integer(int64) :: n, run, lo, mid, hi, left, right, slot

      n = size(keys, 1, kind=int64)
      allocate (order(n), merged(n))
      do lo = 1_int64, n
         order(lo) = lo
      end do

      run = 1_int64
      do while (run < n)
         lo = 1_int64
         do while (lo <= n)
            mid = min(lo + run - 1_int64, n)
            hi = min(lo + 2_int64*run - 1_int64, n)
            left = lo
            right = mid + 1_int64
            do slot = lo, hi
               if (left > mid) then
                  merged(slot) = order(right)
                  right = right + 1_int64
               else if (right > hi) then
                  merged(slot) = order(left)
                  left = left + 1_int64
               else if (key_before(keys(order(right), :), keys(order(left), :))) then
                  merged(slot) = order(right)
                  right = right + 1_int64
               else
                  merged(slot) = order(left)
                  left = left + 1_int64
               end if
            end do
            lo = lo + 2_int64*run
         end do
         order = merged
         run = 2_int64*run
      end do
   end subroutine order_keys

   pure function has_key(keys, order, probe) result(found)
      !! Whether `probe` equals a row of `keys`, which `order` sorts
      integer, intent(in) :: keys(:, :), probe(:)
      integer(int64), intent(in) :: order(:)
      logical :: found

      integer(int64) :: lo, hi, mid

      found = .false.
      lo = 1_int64
      hi = size(order, kind=int64)
      do while (lo <= hi)
         mid = lo + (hi - lo)/2_int64
         if (key_before(keys(order(mid), :), probe)) then
            lo = mid + 1_int64
         else if (key_before(probe, keys(order(mid), :))) then
            hi = mid - 1_int64
         else
            found = .true.
            return
         end if
      end do
   end function has_key

   pure function has_equal_rows(keys, order) result(equal)
      !! Whether two rows of `keys`, sorted by `order`, name the same term
      integer, intent(in) :: keys(:, :)
      integer(int64), intent(in) :: order(:)
      logical :: equal

      integer(int64) :: i

      equal = .false.
      do i = 1_int64, size(order, kind=int64) - 1_int64
         if (.not. key_before(keys(order(i), :), keys(order(i + 1_int64), :))) then
            equal = .true.
            return
         end if
      end do
   end function has_equal_rows

   pure function absolute_key(key) result(probe)
      !! The key with every ghost made real again
      integer, intent(in) :: key(:)
      integer :: probe(size(key))

      integer, allocatable :: entries(:)

      entries = pack(abs(key), key /= 0)
      call sort(entries)
      probe = 0
      probe(1:size(entries)) = entries
   end function absolute_key

   pure function without_entry(key, p) result(probe)
      !! The key with its `p`-th entry removed
      integer, intent(in) :: key(:)
      integer, intent(in) :: p
      integer :: probe(size(key))

      probe = 0
      probe(1:p - 1) = key(1:p - 1)
      probe(p:size(key) - 1) = key(p + 1:size(key))
   end function without_entry

   pure function ghosted_entry(key, p) result(probe)
      !! The key with its `p`-th entry made a ghost
      integer, intent(in) :: key(:)
      integer, intent(in) :: p
      integer :: probe(size(key))

      probe = key
      probe(p) = -key(p)
      call sort(probe(1:count(key /= 0)))
   end function ghosted_entry

   pure function full_basis_key(monomer, n_monomers, width) result(probe)
      !! Monomer `monomer` real and every other monomer of the system ghosted
      integer, intent(in) :: monomer, n_monomers, width
      integer :: probe(width)

      integer :: m

      probe = 0
      probe(1:n_monomers) = [(-m, m=1, n_monomers)]
      probe(monomer) = monomer
      call sort(probe(1:n_monomers))
   end function full_basis_key

   pure function names_whole_system(key, n_monomers) result(whole)
      !! Whether the key names each of the system's monomers exactly once
      integer, intent(in) :: key(:)
      integer, intent(in) :: n_monomers
      logical :: whole

      integer :: m
      integer :: named(size(key))

      named = absolute_key(key)
      whole = all(named(1:n_monomers) == [(m, m=1, n_monomers)])
      if (size(key) > n_monomers) whole = whole .and. all(named(n_monomers + 1:) == 0)
   end function names_whole_system

   pure function fragment_size_of(row) result(n)
      !! How many monomers a polymer row names, padding excluded
      !!
      !! Rows are zero-padded to the widest fragment, so the size is the count
      !! of non-zero entries. Non-zero and not positive: a negative entry is a
      !! monomer present as *ghost centres* -- its atoms and basis functions
      !! without its nucleus -- which the fragment contains and must be sized
      !! for.
      integer(default_int), intent(in) :: row(:)
      integer(default_int) :: n

      n = count(row /= 0)
   end function fragment_size_of

   recursive subroutine generate_fragment_list(monomers, max_level, polymers, count)
      !! Append every combination of 2 to `max_level` monomers to `polymers`
      !!
      !! Monomers are not generated here -- the caller writes those rows first
      !! and passes their number in `count`, which this advances.
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
      !! Generate all combinations of size `r` from the first `n` of `arr`
      integer(default_int), intent(in) :: arr(:)
      integer(default_int), intent(in) :: n, r
      integer(default_int), intent(inout) :: out_array(:, :)
      integer(int64), intent(inout) :: count
      integer(default_int) :: data(r)
      call combine_util(arr, n, r, 1, data, 1, out_array, count)
   end subroutine combine

   recursive subroutine combine_util(arr, n, r, index, data, i, out_array, count)
      !! One level of the combination recursion
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
      !! Print the combinations held in `out_array`, one line each
      integer(default_int), intent(in) ::  max_len
      integer(default_int), intent(in) :: out_array(:, :)
      integer(int64), intent(in) :: count
      integer(int64) :: i
      integer(default_int) :: j

      character(len=MAX_LINE_LENGTH) :: line

      ! Assembled in `line` and emitted once per combination: a log record is a
      ! whole line or it is nothing.
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
      !! Closest approach between different monomers, per fragment, in Angstrom
      !!
      !! Zero for a monomer row.
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
