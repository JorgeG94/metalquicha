!! Occupation-restricted active spaces: the partition, and where a determinant lives
module mqc_ormas_space
   !! A CAS says how many electrons go into how many orbitals and stops there.
   !! ORMAS cuts the active orbitals into consecutive subspaces and puts a
   !! window on the electron count of each, which is enough to express RAS, a
   !! truncated CI, several non-communicating active spaces, or a fragment
   !! model with limited charge transfer -- all as the same object.
   !!
   !! The whole restriction turns out to be one bitmap. Enumerate, per spin, the
   !! *occupation classes* -- the ways the electrons of that spin can be spread
   !! over the subspaces -- and then a determinant is allowed exactly when the
   !! two classes its strings belong to sum, subspace by subspace, into every
   !! window. Nothing depends on which string within a class, so the test is
   !! made once per pair of classes and never again. That is `compatible` below,
   !! and everything else in this module is bookkeeping around it.
   !!
   !! The consequence worth stating plainly: the CI vector stops being a matrix.
   !! Some class pairs are allowed and some are not, so the alpha-by-beta
   !! rectangle acquires holes. What is left is a *list of dense rectangles* --
   !! one per allowed class pair -- laid out alpha-major, and the address of a
   !! determinant is a base for its alpha string plus an offset for its beta
   !! class plus its beta string's rank. A CAS is the case where there is one
   !! class of each spin and the list is a single rectangle, which is why the
   !! formula here reduces to the obvious one and must be checked that it does.
   !!
   !! No point-group symmetry. It is what makes this layout ragged rather than
   !! blocked -- with irreps, the row belonging to an alpha string depends on
   !! that string and not merely on its class -- and this code has no orbital
   !! symmetry anywhere to draw on. The seam for adding it later is narrow and
   !! deliberate: `alpha_base` is stored per string, not recomputed from its
   !! class, and `block_offset` is a table rather than an expression. Both keep
   !! their shape and their consumers when an irrep index appears.
   !!
   !! Like `mqc_determinants`, everything here is checkable without physics.
   !! The class lists are integer compositions, the counts are products of
   !! binomial coefficients, and the address is a bijection onto `1 ..
   !! n_determinants` that has to be dense and has to round trip.
   use pic_types, only: int64, int32
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: MAX_SUBSPACES
   public :: ormas_space_t
   public :: build_ormas_space
   public :: determinant_address
   public :: describes_a_cas

   integer, parameter :: MAX_SUBSPACES = 16
      !! Far above anything meaningful -- the cost of a subspace is another
      !! factor in the class count, and partitions past a handful stop being
      !! interpretable long before they stop being computable. GAMESS allows 50
      !! and checks nothing; a partition that needs more than sixteen is a
      !! mistake worth catching rather than a calculation worth running.

   type :: ormas_space_t
      !! A partition of the active orbitals, and the addressing it implies
      !!
      !! Built once by `build_ormas_space` and read-only afterwards. Every
      !! index is 1-based; `alpha_base` and `block_offset` are offsets and so
      !! count from zero, which is the one place the two conventions meet and
      !! the one place to look first when an address comes out wrong.
      integer :: n_subspaces = 0
      integer :: n_alpha = 0, n_beta = 0
      integer :: n_active_orbitals = 0

      integer, allocatable :: first_orbital(:)
         !! Active orbital each subspace starts at, ascending
      integer, allocatable :: n_orbitals(:)
         !! Orbitals per subspace, derived from `first_orbital`

      integer, allocatable :: min_electrons(:), max_electrons(:)
         !! The windows, **after tightening** -- see `tighten_windows`. These
         !! are the numbers the compatibility test uses, and they can differ
         !! from what the user wrote.
      integer, allocatable :: alpha_min(:), alpha_max(:)
      integer, allocatable :: beta_min(:), beta_max(:)
         !! Per-spin bounds implied by the windows. Necessary but *not*
         !! sufficient: they exist so the two spins can be enumerated
         !! independently, and the real constraint is reimposed by
         !! `compatible`. A class list built from these alone is too generous.

      integer :: n_alpha_classes = 0, n_beta_classes = 0
      integer, allocatable :: alpha_class(:, :), beta_class(:, :)
         !! `(subspace, class)` -- the occupation vectors, in enumeration order
      integer(int64), allocatable :: alpha_class_size(:), beta_class_size(:)
         !! Strings in each class: the product over subspaces of `C(orbitals,
         !! electrons)`, since within a class the subspaces are independent
      integer(int64), allocatable :: alpha_offset(:), beta_offset(:)
         !! Running string totals, size `n_classes + 1`. Class `g` owns global
         !! string indices `offset(g) + 1 .. offset(g + 1)`.
      integer, allocatable :: alpha_string_class(:), beta_string_class(:)
         !! Which class each string belongs to. Stored rather than searched:
         !! it is read in the innermost loop of everything downstream.

      logical, allocatable :: compatible(:, :)
         !! `(beta class, alpha class)` -- the entire ORMAS restriction. Beta
         !! first because that is the order it is consumed in.

      integer(int64), allocatable :: row_length(:)
         !! Determinants carried by one alpha string of each class
      integer(int64), allocatable :: class_base(:)
         !! Determinants preceding all of an alpha class
      integer(int64), allocatable :: block_offset(:, :)
         !! `(beta class, alpha class)` -- where that beta class starts inside
         !! an alpha string's row. Meaningless where `compatible` is false, and
         !! left at zero there.
      integer(int64), allocatable :: alpha_base(:)
         !! Determinants preceding every determinant of each alpha *string*

      integer(int64) :: n_determinants = 0_int64
   contains
      procedure :: destroy => ormas_space_destroy
   end type ormas_space_t

contains

   pure function binomial(n, k) result(c)
      !! `C(n, k)`, exactly, in 64-bit integers
      !!
      !! One factor at a time, multiplying before dividing, so the running
      !! value is always the binomial coefficient of a smaller problem and
      !! never exceeds the answer.
      integer, intent(in) :: n, k
      integer(int64) :: c

      integer :: i, kk

      c = 0_int64
      if (k < 0 .or. k > n .or. n < 0) return
      kk = min(k, n - k)
      c = 1_int64
      do i = 1, kk
         c = c*int(n - kk + i, int64)/int(i, int64)
      end do
   end function binomial

   pure function class_size(n_orbitals, occupation) result(count)
      !! Strings realising one occupation vector
      !!
      !! The subspaces are independent given the class, so this is a product.
      integer, intent(in) :: n_orbitals(:), occupation(:)
      integer(int64) :: count

      integer :: k

      count = 1_int64
      do k = 1, size(occupation)
         count = count*binomial(n_orbitals(k), occupation(k))
      end do
   end function class_size

   pure subroutine first_composition(capacity, minimum, total, occupation)
      !! The first way to spread `total` particles over the boxes
      !!
      !! Every box starts at its minimum, then the surplus is poured in from
      !! the left until it runs out. Filling greedily from the left is what
      !! makes the enumeration order descending-lexicographic, which
      !! `next_composition` then walks and which the addressing assumes.
      integer, intent(in) :: capacity(:), minimum(:)
      integer, intent(in) :: total
      integer, intent(out) :: occupation(:)

      integer :: k, surplus, room

      occupation = minimum
      surplus = total - sum(minimum)
      do k = 1, size(capacity)
         room = min(surplus, capacity(k) - minimum(k))
         occupation(k) = minimum(k) + room
         surplus = surplus - room
      end do
   end subroutine first_composition

   pure subroutine next_composition(capacity, minimum, total, occupation, exhausted)
      !! The next way, or a flag saying there is none
      !!
      !! Find the rightmost box that can give a particle to some box on its
      !! right, move one out of it, and repack everything to its right greedily
      !! from the left. That is the successor in descending-lexicographic
      !! order, and it visits every composition exactly once.
      integer, intent(in) :: capacity(:), minimum(:)
      integer, intent(in) :: total
      integer, intent(inout) :: occupation(:)
      logical, intent(out) :: exhausted

      integer :: n_boxes, k, m, remaining, room

      n_boxes = size(capacity)
      exhausted = .true.

      do k = n_boxes - 1, 1, -1
         if (occupation(k) <= minimum(k)) cycle
         if (.not. any(occupation(k + 1:) < capacity(k + 1:))) cycle

         occupation(k) = occupation(k) - 1
         remaining = total - sum(occupation(1:k))
         do m = k + 1, n_boxes
            room = min(remaining - sum(minimum(m + 1:)), capacity(m))
            occupation(m) = max(minimum(m), room)
            remaining = remaining - occupation(m)
         end do
         exhausted = .false.
         return
      end do
   end subroutine next_composition

   subroutine enumerate_classes(capacity, minimum, total, classes, n_classes)
      !! Every occupation vector allowed by the per-spin bounds
      !!
      !! Counted first and then filled, because the count is wanted on its own
      !! and running the walk twice is cheaper than growing an array.
      integer, intent(in) :: capacity(:), minimum(:)
      integer, intent(in) :: total
      integer, allocatable, intent(out) :: classes(:, :)
      integer, intent(out) :: n_classes

      integer :: occupation(size(capacity))
      integer :: n_boxes
      logical :: exhausted

      n_boxes = size(capacity)
      n_classes = 0

      if (total < sum(minimum) .or. total > sum(capacity)) then
         allocate (classes(n_boxes, 0))
         return
      end if

      call first_composition(capacity, minimum, total, occupation)
      do
         n_classes = n_classes + 1
         call next_composition(capacity, minimum, total, occupation, exhausted)
         if (exhausted) exit
      end do

      allocate (classes(n_boxes, n_classes))
      call first_composition(capacity, minimum, total, occupation)
      n_classes = 0
      do
         n_classes = n_classes + 1
         classes(:, n_classes) = occupation
         call next_composition(capacity, minimum, total, occupation, exhausted)
         if (exhausted) exit
      end do
   end subroutine enumerate_classes

   subroutine validate_partition(space, n_electrons, error)
      !! The ways a partition can fail to describe anything
      type(ormas_space_t), intent(in) :: space
      integer, intent(in) :: n_electrons
      type(error_t), intent(inout) :: error

      integer :: k

      if (error%has_error()) return

      if (space%n_subspaces < 1 .or. space%n_subspaces > MAX_SUBSPACES) then
         call error%set(ERROR_VALIDATION, "an active space cut into "// &
                        to_char(space%n_subspaces)//" subspaces is not one this "// &
                        "code will build; the limit is "//to_char(MAX_SUBSPACES)//".")
         return
      end if

      do k = 1, space%n_subspaces
         if (space%n_orbitals(k) < 1) then
            call error%set(ERROR_VALIDATION, "subspace "//to_char(k)//" has no "// &
                           "orbitals in it. The subspaces are consecutive ranges of "// &
                           "the active orbitals, so their starting orbitals must "// &
                           "ascend and the last must leave room.")
            return
         end if
         if (space%min_electrons(k) > space%max_electrons(k)) then
            call error%set(ERROR_VALIDATION, "subspace "//to_char(k)//" is asked to "// &
                           "hold at least "//to_char(space%min_electrons(k))// &
                           " electrons and at most "//to_char(space%max_electrons(k))// &
                           ", which no subspace can do.")
            return
         end if
         if (space%max_electrons(k) > 2*space%n_orbitals(k)) then
            call error%set(ERROR_VALIDATION, "subspace "//to_char(k)//" is allowed "// &
                           to_char(space%max_electrons(k))//" electrons but has only "// &
                           to_char(space%n_orbitals(k))//" orbitals, which hold "// &
                           to_char(2*space%n_orbitals(k))//".")
            return
         end if
      end do

      if (sum(space%min_electrons) > n_electrons) then
         call error%set(ERROR_VALIDATION, "the subspace minima add to "// &
                        to_char(sum(space%min_electrons))//" electrons and there are "// &
                        to_char(n_electrons)//" to place, so no determinant satisfies "// &
                        "all of them at once.")
         return
      end if
      if (sum(space%max_electrons) < n_electrons) then
         call error%set(ERROR_VALIDATION, "the subspace maxima add to "// &
                        to_char(sum(space%max_electrons))//" electrons and there are "// &
                        to_char(n_electrons)//" to place, so they do not all fit.")
         return
      end if
   end subroutine validate_partition

   pure subroutine tighten_windows(space, n_electrons)
      !! Shrink the windows to what the other subspaces actually leave available
      !!
      !! A window can promise more room than the rest of the partition can
      !! spare: if every other subspace must hold at least so many electrons,
      !! this one cannot exceed what is left over however large its own maximum
      !! reads. Tightening first means the compatibility test sees the real
      !! bounds, and it is why the numbers reported back may not be the numbers
      !! that were asked for.
      !!
      !! Order matters -- the minima are tightened against the *already*
      !! tightened maxima.
      type(ormas_space_t), intent(inout) :: space
      integer, intent(in) :: n_electrons

      integer :: k, spare, reachable

      spare = n_electrons - sum(space%min_electrons)
      do k = 1, space%n_subspaces
         space%max_electrons(k) = space%min_electrons(k) + &
                                  min(space%max_electrons(k) - space%min_electrons(k), spare)
      end do

      do k = 1, space%n_subspaces
         reachable = n_electrons - (sum(space%max_electrons) - space%max_electrons(k))
         space%min_electrons(k) = max(space%min_electrons(k), reachable)
      end do
   end subroutine tighten_windows

   pure subroutine derive_spin_bounds(space)
      !! The loosest per-spin bounds the windows imply
      !!
      !! The windows constrain alpha and beta together, which would force the
      !! two spins to be enumerated jointly. These bounds are the tightest
      !! box that contains the joint constraint and factorises, so the strings
      !! can be generated one spin at a time; what they let through that the
      !! real constraint does not is caught later, once per class pair, by
      !! `compatible`. Using them *as* the constraint would be wrong.
      type(ormas_space_t), intent(inout) :: space

      integer :: k, room, forced_by_window, forced_by_others

      do k = 1, space%n_subspaces
         room = min(space%n_orbitals(k), space%max_electrons(k))
         space%alpha_max(k) = min(room, space%n_alpha)
         space%beta_max(k) = min(room, space%n_beta)
      end do

      do k = 1, space%n_subspaces
         forced_by_window = max(0, space%min_electrons(k) - space%beta_max(k))
         forced_by_others = space%n_alpha - (sum(space%alpha_max) - space%alpha_max(k))
         space%alpha_min(k) = max(forced_by_window, forced_by_others)

         forced_by_window = max(0, space%min_electrons(k) - space%alpha_max(k))
         forced_by_others = space%n_beta - (sum(space%beta_max) - space%beta_max(k))
         space%beta_min(k) = max(forced_by_window, forced_by_others)
      end do
   end subroutine derive_spin_bounds

   subroutine class_string_tables(classes, n_orbitals, sizes, offsets, string_class)
      !! Per-class string counts, their running totals, and the inverse map
      integer, intent(in) :: classes(:, :), n_orbitals(:)
      integer(int64), allocatable, intent(out) :: sizes(:), offsets(:)
      integer, allocatable, intent(out) :: string_class(:)

      integer :: n_classes, g
      integer(int64) :: s

      n_classes = size(classes, 2)
      allocate (sizes(n_classes), offsets(n_classes + 1))

      offsets(1) = 0_int64
      do g = 1, n_classes
         sizes(g) = class_size(n_orbitals, classes(:, g))
         offsets(g + 1) = offsets(g) + sizes(g)
      end do

      allocate (string_class(offsets(n_classes + 1)))
      do g = 1, n_classes
         do s = offsets(g) + 1, offsets(g + 1)
            string_class(s) = g
         end do
      end do
   end subroutine class_string_tables

   subroutine build_addressing(space, error)
      !! Turn the compatibility bitmap into an address for every determinant
      !!
      !! Alpha-major: each alpha string owns a contiguous run of determinants,
      !! and inside that run the beta classes appear in order, each contributing
      !! a dense block of its own strings. So three tables suffice -- how long
      !! a row is, where each beta class starts within it, and where each alpha
      !! string's row starts overall.
      type(ormas_space_t), intent(inout) :: space
      type(error_t), intent(inout) :: error

      integer :: ga, gb, a
      integer(int64) :: running, total

      allocate (space%row_length(space%n_alpha_classes))
      allocate (space%class_base(space%n_alpha_classes))
      allocate (space%block_offset(space%n_beta_classes, space%n_alpha_classes))

      space%block_offset = 0_int64

      do ga = 1, space%n_alpha_classes
         running = 0_int64
         do gb = 1, space%n_beta_classes
            if (.not. space%compatible(gb, ga)) cycle
            space%block_offset(gb, ga) = running
            running = running + space%beta_class_size(gb)
         end do
         space%row_length(ga) = running
      end do

      total = 0_int64
      do ga = 1, space%n_alpha_classes
         space%class_base(ga) = total
         total = total + space%alpha_class_size(ga)*space%row_length(ga)
      end do
      space%n_determinants = total

      if (total > int(huge(1_int32), int64)) then
         call error%set(ERROR_VALIDATION, "this partition has "//to_char(total)// &
                        " determinants, past what a default integer can address and "// &
                        "far past what the expansion could be stored for. Tighten the "// &
                        "subspace windows or shrink the active space.")
         return
      end if

      allocate (space%alpha_base(space%alpha_offset(space%n_alpha_classes + 1)))
      do a = 1, int(space%alpha_offset(space%n_alpha_classes + 1))
         ga = space%alpha_string_class(a)
         space%alpha_base(a) = space%class_base(ga) + &
                               (int(a, int64) - space%alpha_offset(ga) - 1_int64)* &
                               space%row_length(ga)
      end do
   end subroutine build_addressing

   subroutine build_ormas_space(first_orbital, n_active_orbitals, n_alpha, n_beta, &
                                min_electrons, max_electrons, space, error)
      !! Assemble a partition from what a deck can state
      !!
      !! `first_orbital` gives the active orbital each subspace begins at, so
      !! its length is the number of subspaces and its first entry is 1. The
      !! windows are on the *total* electron count of a subspace, both spins
      !! together, which is what makes them independent of how the spins are
      !! distributed and what lets the whole restriction reduce to a bitmap.
      integer, intent(in) :: first_orbital(:)
      integer, intent(in) :: n_active_orbitals, n_alpha, n_beta
      integer, intent(in) :: min_electrons(:), max_electrons(:)
      type(ormas_space_t), intent(out) :: space
      type(error_t), intent(inout) :: error

      integer :: n_subspaces, k, ga, gb
      logical :: allowed

      if (error%has_error()) return

      n_subspaces = size(first_orbital)
      if (size(min_electrons) /= n_subspaces .or. size(max_electrons) /= n_subspaces) then
         call error%set(ERROR_VALIDATION, "the partition names "//to_char(n_subspaces)// &
                        " subspaces but gives "//to_char(size(min_electrons))// &
                        " minima and "//to_char(size(max_electrons))//" maxima.")
         return
      end if
      if (n_subspaces < 1 .or. n_subspaces > MAX_SUBSPACES) then
         call error%set(ERROR_VALIDATION, "an active space cut into "// &
                        to_char(n_subspaces)//" subspaces is not one this code will "// &
                        "build; the limit is "//to_char(MAX_SUBSPACES)//".")
         return
      end if

      space%n_subspaces = n_subspaces
      space%n_alpha = n_alpha
      space%n_beta = n_beta
      space%n_active_orbitals = n_active_orbitals

      allocate (space%first_orbital(n_subspaces), space%n_orbitals(n_subspaces))
      allocate (space%min_electrons(n_subspaces), space%max_electrons(n_subspaces))
      allocate (space%alpha_min(n_subspaces), space%alpha_max(n_subspaces))
      allocate (space%beta_min(n_subspaces), space%beta_max(n_subspaces))

      space%first_orbital = first_orbital
      space%min_electrons = min_electrons
      space%max_electrons = max_electrons

      do k = 1, n_subspaces - 1
         space%n_orbitals(k) = first_orbital(k + 1) - first_orbital(k)
      end do
      space%n_orbitals(n_subspaces) = n_active_orbitals - first_orbital(n_subspaces) + 1

      if (first_orbital(1) /= 1) then
         call error%set(ERROR_VALIDATION, "the first subspace starts at active orbital "// &
                        to_char(first_orbital(1))//"; the subspaces have to cover the "// &
                        "active orbitals, so it starts at 1.")
         return
      end if

      call validate_partition(space, n_alpha + n_beta, error)
      if (error%has_error()) return

      call tighten_windows(space, n_alpha + n_beta)
      call derive_spin_bounds(space)

      call enumerate_classes(space%alpha_max, space%alpha_min, n_alpha, &
                             space%alpha_class, space%n_alpha_classes)
      call enumerate_classes(space%beta_max, space%beta_min, n_beta, &
                             space%beta_class, space%n_beta_classes)

      if (space%n_alpha_classes == 0 .or. space%n_beta_classes == 0) then
         call error%set(ERROR_VALIDATION, "no way of spreading "//to_char(n_alpha)// &
                        " alpha and "//to_char(n_beta)//" beta electrons over these "// &
                        "subspaces satisfies their windows.")
         return
      end if

      call class_string_tables(space%alpha_class, space%n_orbitals, &
                               space%alpha_class_size, space%alpha_offset, &
                               space%alpha_string_class)
      call class_string_tables(space%beta_class, space%n_orbitals, &
                               space%beta_class_size, space%beta_offset, &
                               space%beta_string_class)

      allocate (space%compatible(space%n_beta_classes, space%n_alpha_classes))
      do ga = 1, space%n_alpha_classes
         do gb = 1, space%n_beta_classes
            allowed = .true.
            do k = 1, n_subspaces
               associate (occupancy => space%alpha_class(k, ga) + space%beta_class(k, gb))
                  if (occupancy < space%min_electrons(k)) allowed = .false.
                  if (occupancy > space%max_electrons(k)) allowed = .false.
               end associate
            end do
            space%compatible(gb, ga) = allowed
         end do
      end do

      if (.not. any(space%compatible)) then
         call error%set(ERROR_VALIDATION, "no pairing of an alpha and a beta "// &
                        "occupation satisfies the subspace windows, so the space is "// &
                        "empty. The windows are on both spins together.")
         return
      end if

      call build_addressing(space, error)
   end subroutine build_ormas_space

   pure function determinant_address(space, alpha_string, beta_string) result(address)
      !! Where the determinant made of these two strings sits
      !!
      !! The three terms are the alpha string's row, the beta class's block
      !! within that row, and the beta string's rank within its class. Valid
      !! only for a pair whose classes are compatible -- an incompatible pair
      !! names no determinant, and this returns a plausible-looking number for
      !! it rather than a signal, which is why the caller gates on
      !! `compatible` and not on the answer.
      type(ormas_space_t), intent(in) :: space
      integer, intent(in) :: alpha_string, beta_string
      integer(int64) :: address

      integer :: ga, gb

      ga = space%alpha_string_class(alpha_string)
      gb = space%beta_string_class(beta_string)

      address = space%alpha_base(alpha_string) &
                + space%block_offset(gb, ga) &
                + (int(beta_string, int64) - space%beta_offset(gb))
   end function determinant_address

   pure function describes_a_cas(space) result(is_cas)
      !! Whether this partition restricts anything at all
      !!
      !! Windows that happen to allow every distribution leave a complete
      !! active space, however many subspaces are drawn on it. The test is
      !! whether the determinant count matches the unrestricted one, not
      !! whether there is more than one subspace -- the distinction matters
      !! downstream, where a genuine CAS makes rotations between active
      !! orbitals redundant and an ORMAS space does not.
      type(ormas_space_t), intent(in) :: space
      logical :: is_cas

      is_cas = space%n_determinants == &
         binomial(space%n_active_orbitals, space%n_alpha)* &
         binomial(space%n_active_orbitals, space%n_beta)
   end function describes_a_cas

   pure subroutine ormas_space_destroy(this)
      !! Give back everything the partition holds
      class(ormas_space_t), intent(inout) :: this

      if (allocated(this%first_orbital)) deallocate (this%first_orbital)
      if (allocated(this%n_orbitals)) deallocate (this%n_orbitals)
      if (allocated(this%min_electrons)) deallocate (this%min_electrons)
      if (allocated(this%max_electrons)) deallocate (this%max_electrons)
      if (allocated(this%alpha_min)) deallocate (this%alpha_min)
      if (allocated(this%alpha_max)) deallocate (this%alpha_max)
      if (allocated(this%beta_min)) deallocate (this%beta_min)
      if (allocated(this%beta_max)) deallocate (this%beta_max)
      if (allocated(this%alpha_class)) deallocate (this%alpha_class)
      if (allocated(this%beta_class)) deallocate (this%beta_class)
      if (allocated(this%alpha_class_size)) deallocate (this%alpha_class_size)
      if (allocated(this%beta_class_size)) deallocate (this%beta_class_size)
      if (allocated(this%alpha_offset)) deallocate (this%alpha_offset)
      if (allocated(this%beta_offset)) deallocate (this%beta_offset)
      if (allocated(this%alpha_string_class)) deallocate (this%alpha_string_class)
      if (allocated(this%beta_string_class)) deallocate (this%beta_string_class)
      if (allocated(this%compatible)) deallocate (this%compatible)
      if (allocated(this%row_length)) deallocate (this%row_length)
      if (allocated(this%class_base)) deallocate (this%class_base)
      if (allocated(this%block_offset)) deallocate (this%block_offset)
      if (allocated(this%alpha_base)) deallocate (this%alpha_base)

      this%n_subspaces = 0
      this%n_alpha_classes = 0
      this%n_beta_classes = 0
      this%n_determinants = 0_int64
   end subroutine ormas_space_destroy

end module mqc_ormas_space
