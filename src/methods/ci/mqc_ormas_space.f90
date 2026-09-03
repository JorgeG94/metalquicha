!! Occupation-restricted active spaces: the partition, and where a determinant lives
module mqc_ormas_space
   !! ORMAS cuts the active orbitals into consecutive subspaces and puts a
   !! window on the electron count of each, which expresses RAS, a truncated CI,
   !! several non-communicating active spaces, or a fragment model with limited
   !! charge transfer as one object.
   !!
   !! The whole restriction is one bitmap. Per spin, the *occupation classes*
   !! are the ways that spin's electrons can be spread over the subspaces; a
   !! determinant is allowed exactly when the two classes its strings belong to
   !! sum, subspace by subspace, into every window. That is `compatible` below.
   !!
   !! So the CI vector is not a matrix but a *list of dense rectangles*, one per
   !! allowed class pair, laid out alpha-major: the address of a determinant is
   !! a base for its alpha string, plus an offset for its beta class, plus its
   !! beta string's rank. A CAS is the single-rectangle case.
   !!
   !! No point-group symmetry: nothing here is blocked by irrep, and
   !! `alpha_base` is stored per string rather than recomputed from its class.
   use pic_types, only: int64, int32
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_determinants, only: generate_strings, string_address
   use mqc_math_utils, only: binomial
   implicit none
   private

   public :: MAX_SUBSPACES
   public :: ormas_space_t
   public :: build_ormas_space
   public :: determinant_address
   public :: describes_a_cas
   public :: ormas_strings
   public :: ormas_string_address
   public :: ormas_closure_t
   public :: build_ormas_closure
   public :: closure_address

   integer, parameter :: MAX_SUBSPACES = 16
      !! Each subspace is another factor in the class count, so a partition
      !! needing more than this is a mistake rather than a calculation.

   type :: ormas_space_t
      !! A partition of the active orbitals, and the addressing it implies
      !!
      !! Built once by `build_ormas_space` and read-only afterwards. Every
      !! index is 1-based; `alpha_base` and `block_offset` are offsets and so
      !! count from zero.
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
      integer :: n_alpha_classes = 0, n_beta_classes = 0
      integer, allocatable :: alpha_class(:, :), beta_class(:, :)
         !! `(subspace, class)` -- every way that spin's electrons can be spread
         !! over the subspaces, in enumeration order, whether or not the windows
         !! leave any use for it.
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
         !! `(beta class, alpha class)` -- the entire ORMAS restriction

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

   type :: ormas_closure_t
      !! The space plus one excitation, and where a determinant sits in it
      !!
      !! An excitation may take a determinant out of the space and a second
      !! bring it back, so the intermediate vector of a sigma build lives on the
      !! space widened by one excitation, not on the space. Truncating it to the
      !! space loses those terms silently.
      !!
      !! Same class lists, same strings and the same three addressing tables as
      !! `ormas_space_t`; only which pairs are in it differs.
      logical, allocatable :: present(:, :)
         !! `(beta class, alpha class)` -- reachable in one excitation, or in
         !! the space already
      integer(int64), allocatable :: row_length(:), class_base(:)
      integer(int64), allocatable :: block_offset(:, :), alpha_base(:)
      integer(int64) :: n_determinants = 0_int64
   contains
      procedure :: destroy => ormas_closure_destroy
   end type ormas_closure_t

contains

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
      !! Boxes start at their minima and the surplus fills from the left, which
      !! puts the enumeration in descending-lexicographic order -- what
      !! `next_composition` walks and what the addressing assumes.
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
      !! The successor in descending-lexicographic order; the walk visits every
      !! composition exactly once.
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

      ! TODO(mqc): `build_ormas_space` makes this same check, with the same
      ! message, before it calls here, so this copy can never fire.
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
      !! spare, so the bounds stored afterwards may not be the ones that were
      !! asked for. The minima are tightened against the *already* tightened
      !! maxima, so the order matters.
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

   subroutine layout_addressing(grid, alpha_size, beta_size, alpha_offset, &
                                alpha_string_class, row_length, class_base, &
                                block_offset, alpha_base, n_determinants, error)
      !! Turn a grid of allowed class pairs into an address for every pair
      !!
      !! Alpha-major: each alpha string owns a contiguous run, and inside it the
      !! beta classes appear in order, each contributing a dense block of its own
      !! strings. Three tables say it all -- how long a row is, where each beta
      !! class starts within one, and where each alpha string's row begins.
      logical, intent(in) :: grid(:, :)
      integer(int64), intent(in) :: alpha_size(:), beta_size(:), alpha_offset(:)
      integer, intent(in) :: alpha_string_class(:)
      integer(int64), allocatable, intent(out) :: row_length(:), class_base(:)
      integer(int64), allocatable, intent(out) :: block_offset(:, :), alpha_base(:)
      integer(int64), intent(out) :: n_determinants
      type(error_t), intent(inout) :: error

      integer :: n_alpha_classes, n_beta_classes, ga, gb, a
      integer(int64) :: running, total

      n_alpha_classes = size(alpha_size)
      n_beta_classes = size(beta_size)

      allocate (row_length(n_alpha_classes), class_base(n_alpha_classes))
      allocate (block_offset(n_beta_classes, n_alpha_classes))
      block_offset = 0_int64

      do ga = 1, n_alpha_classes
         running = 0_int64
         do gb = 1, n_beta_classes
            if (.not. grid(gb, ga)) cycle
            block_offset(gb, ga) = running
            running = running + beta_size(gb)
         end do
         row_length(ga) = running
      end do

      total = 0_int64
      do ga = 1, n_alpha_classes
         class_base(ga) = total
         total = total + alpha_size(ga)*row_length(ga)
      end do
      n_determinants = total

      if (total > int(huge(1_int32), int64)) then
         call error%set(ERROR_VALIDATION, "this partition has "//to_char(total)// &
                        " determinants, past what a default integer can address and "// &
                        "far past what the expansion could be stored for. Tighten the "// &
                        "subspace windows or shrink the active space.")
         return
      end if

      allocate (alpha_base(size(alpha_string_class)))
      do a = 1, size(alpha_string_class)
         ga = alpha_string_class(a)
         alpha_base(a) = class_base(ga) + &
                         (int(a, int64) - alpha_offset(ga) - 1_int64)*row_length(ga)
      end do
   end subroutine layout_addressing

   subroutine build_grid(space)
      !! Which pairs of occupation classes the windows allow
      type(ormas_space_t), intent(inout) :: space

      integer :: ga, gb, k
      logical :: allowed

      if (allocated(space%compatible)) deallocate (space%compatible)
      allocate (space%compatible(space%n_beta_classes, space%n_alpha_classes))

      do ga = 1, space%n_alpha_classes
         do gb = 1, space%n_beta_classes
            allowed = .true.
            do k = 1, space%n_subspaces
               associate (occupancy => space%alpha_class(k, ga) + space%beta_class(k, gb))
                  if (occupancy < space%min_electrons(k)) allowed = .false.
                  if (occupancy > space%max_electrons(k)) allowed = .false.
               end associate
            end do
            space%compatible(gb, ga) = allowed
         end do
      end do
   end subroutine build_grid

   subroutine drop_unreachable_classes(space)
      !! Keep the classes the space uses and the ones an excitation can reach
      !!
      !! A class earns its place by appearing in some allowed pair, or by being
      !! one excitation away from one that does. Anything further can hold no
      !! determinant of the space and can be reached by nothing the sigma build
      !! forms.
      !!
      !! **One move, not none.** Pruning to the classes that appear in the space
      !! is wrong here: the sigma intermediate is a vector indexed by
      !! determinant, so every class the closure touches has to exist.
      type(ormas_space_t), intent(inout) :: space

      logical, allocatable :: keep_alpha(:), keep_beta(:)
      integer :: ga, gb

      allocate (keep_alpha(space%n_alpha_classes), keep_beta(space%n_beta_classes))
      keep_alpha = .false.
      keep_beta = .false.

      do ga = 1, space%n_alpha_classes
         do gb = 1, space%n_beta_classes
            if (.not. space%compatible(gb, ga)) cycle
            keep_alpha(ga) = .true.
            keep_beta(gb) = .true.
         end do
      end do

      call widen_by_one_move(space%alpha_class, keep_alpha)
      call widen_by_one_move(space%beta_class, keep_beta)

      call keep_classes(space%alpha_class, keep_alpha, space%n_alpha_classes)
      call keep_classes(space%beta_class, keep_beta, space%n_beta_classes)

      deallocate (keep_alpha, keep_beta)
   end subroutine drop_unreachable_classes

   pure subroutine widen_by_one_move(classes, keep)
      !! Add every class one excitation away from one already kept
      integer, intent(in) :: classes(:, :)
      logical, intent(inout) :: keep(:)

      logical :: seed(size(keep))
      integer :: g, s

      seed = keep
      do g = 1, size(keep)
         if (keep(g)) cycle
         do s = 1, size(keep)
            if (.not. seed(s)) cycle
            if (one_move_apart(classes(:, g), classes(:, s))) then
               keep(g) = .true.
               exit
            end if
         end do
      end do
   end subroutine widen_by_one_move

   subroutine keep_classes(classes, keep, n_classes)
      !! Compact a class list down to the kept entries, in the same order
      integer, allocatable, intent(inout) :: classes(:, :)
      logical, intent(in) :: keep(:)
      integer, intent(out) :: n_classes

      integer, allocatable :: kept(:, :)
      integer :: g

      allocate (kept(size(classes, 1), count(keep)))
      n_classes = 0
      do g = 1, size(keep)
         if (.not. keep(g)) cycle
         n_classes = n_classes + 1
         kept(:, n_classes) = classes(:, g)
      end do
      call move_alloc(kept, classes)
   end subroutine keep_classes

   subroutine build_addressing(space, error)
      !! The space's own addressing, from its compatibility grid
      type(ormas_space_t), intent(inout) :: space
      type(error_t), intent(inout) :: error

      call layout_addressing(space%compatible, space%alpha_class_size, &
                             space%beta_class_size, space%alpha_offset, &
                             space%alpha_string_class, space%row_length, &
                             space%class_base, space%block_offset, space%alpha_base, &
                             space%n_determinants, error)
   end subroutine build_addressing

   subroutine build_ormas_space(first_orbital, n_active_orbitals, n_alpha, n_beta, &
                                min_electrons, max_electrons, space, error)
      !! Assemble a partition from what a deck can state
      !!
      !! `first_orbital` gives the active orbital each subspace begins at, so
      !! its length is the number of subspaces and its first entry is 1. The
      !! windows are on the *total* electron count of a subspace, both spins
      !! together.
      integer, intent(in) :: first_orbital(:)
      integer, intent(in) :: n_active_orbitals, n_alpha, n_beta
      integer, intent(in) :: min_electrons(:), max_electrons(:)
      type(ormas_space_t), intent(out) :: space
      type(error_t), intent(inout) :: error

      integer :: n_subspaces, k

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

      ! Capacity is what a subspace can physically hold of one spin, and there
      ! is no lower bound: the windows constrain the two spins together, so
      ! neither alone is restricted by them. `compatible` decides everything.
      call enumerate_classes(space%n_orbitals, spread(0, 1, n_subspaces), n_alpha, &
                             space%alpha_class, space%n_alpha_classes)
      call enumerate_classes(space%n_orbitals, spread(0, 1, n_subspaces), n_beta, &
                             space%beta_class, space%n_beta_classes)

      if (space%n_alpha_classes == 0 .or. space%n_beta_classes == 0) then
         call error%set(ERROR_VALIDATION, "no way of spreading "//to_char(n_alpha)// &
                        " alpha and "//to_char(n_beta)//" beta electrons over these "// &
                        "subspaces satisfies their windows.")
         return
      end if

      call build_grid(space)
      call drop_unreachable_classes(space)
      call build_grid(space)

      call class_string_tables(space%alpha_class, space%n_orbitals, &
                               space%alpha_class_size, space%alpha_offset, &
                               space%alpha_string_class)
      call class_string_tables(space%beta_class, space%n_orbitals, &
                               space%beta_class_size, space%beta_offset, &
                               space%beta_string_class)

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
      !! only for a compatible pair of classes: an incompatible pair names no
      !! determinant and gets a plausible-looking number rather than a signal,
      !! so the caller checks `compatible` and not the answer.
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
      !! Windows that happen to allow every distribution leave a complete active
      !! space, however many subspaces are drawn on it, so the test is on the
      !! determinant count and not on the number of subspaces. A genuine CAS
      !! makes rotations between active orbitals redundant; an ORMAS space does
      !! not.
      type(ormas_space_t), intent(in) :: space
      logical :: is_cas

      is_cas = space%n_determinants == &
         binomial(space%n_active_orbitals, space%n_alpha)* &
         binomial(space%n_active_orbitals, space%n_beta)
   end function describes_a_cas

   pure function subspace_string(space, subspace, string) result(part)
      !! The bits of a string belonging to one subspace, shifted down to zero
      type(ormas_space_t), intent(in) :: space
      integer, intent(in) :: subspace
      integer(int64), intent(in) :: string
      integer(int64) :: part

      integer :: shift, width

      shift = space%first_orbital(subspace) - 1
      width = space%n_orbitals(subspace)
      part = iand(ishft(string, -shift), ishft(1_int64, width) - 1_int64)
   end function subspace_string

   pure function occupation_of(space, string) result(occupation)
      !! How many electrons a string puts in each subspace
      type(ormas_space_t), intent(in) :: space
      integer(int64), intent(in) :: string
      integer :: occupation(space%n_subspaces)

      integer :: k

      do k = 1, space%n_subspaces
         occupation(k) = popcnt(subspace_string(space, k, string))
      end do
   end function occupation_of

   pure function class_holding(classes, occupation) result(which)
      !! Which class an occupation vector is, or zero if the bounds exclude it
      integer, intent(in) :: classes(:, :), occupation(:)
      integer :: which

      integer :: g

      which = 0
      do g = 1, size(classes, 2)
         if (all(classes(:, g) == occupation)) then
            which = g
            return
         end if
      end do
   end function class_holding

   subroutine strings_of_one_spin(space, classes, offsets, strings, error)
      !! Every string of one spin, in global index order
      !!
      !! Within a class the subspaces are independent, so the strings are the
      !! product of each subspace's own strings, with the *last* subspace
      !! varying fastest -- the order the mixed-radix address in
      !! `within_class_rank` inverts. Inside a subspace the strings come from
      !! `mqc_determinants`, in ascending bit-pattern order.
      type(ormas_space_t), intent(in) :: space
      integer, intent(in) :: classes(:, :)
      integer(int64), intent(in) :: offsets(:)
      integer(int64), allocatable, intent(out) :: strings(:)
      type(error_t), intent(inout) :: error

      integer(int64), allocatable :: part(:)
      integer(int64) :: index, stride, block
      integer :: n_classes, g, k, shift
      integer :: sizes(space%n_subspaces)

      if (error%has_error()) return

      n_classes = size(classes, 2)
      allocate (strings(offsets(n_classes + 1)))
      strings = 0_int64

      do g = 1, n_classes
         do k = 1, space%n_subspaces
            sizes(k) = int(binomial(space%n_orbitals(k), classes(k, g)))
         end do

         stride = 1_int64
         do k = space%n_subspaces, 1, -1
            call generate_strings(space%n_orbitals(k), classes(k, g), part, error)
            if (error%has_error()) return

            shift = space%first_orbital(k) - 1
            ! Repeat this subspace's strings every `stride` determinants, so
            ! subspace `n_subspaces` changes on every step and subspace 1 only
            ! after a full sweep of all the others.
            do index = offsets(g) + 1_int64, offsets(g + 1)
               block = mod((index - offsets(g) - 1_int64)/stride, int(sizes(k), int64))
               strings(index) = ior(strings(index), ishft(part(block + 1_int64), shift))
            end do

            stride = stride*int(sizes(k), int64)
            deallocate (part)
         end do
      end do
   end subroutine strings_of_one_spin

   pure function within_class_rank(space, classes, which, string) result(rank)
      !! Where a string sits among those of its own class, 1-based
      !!
      !! Mixed radix over the subspaces with the last as the least significant
      !! digit, each digit being that subspace's own address.
      type(ormas_space_t), intent(in) :: space
      integer, intent(in) :: classes(:, :)
      integer, intent(in) :: which
      integer(int64), intent(in) :: string
      integer(int64) :: rank

      integer :: k, j
      integer(int64) :: digit, weight

      rank = 1_int64
      do k = 1, space%n_subspaces
         digit = int(string_address(space%n_orbitals(k), classes(k, which), &
                                    subspace_string(space, k, string)), int64) - 1_int64
         weight = 1_int64
         do j = k + 1, space%n_subspaces
            weight = weight*binomial(space%n_orbitals(j), classes(j, which))
         end do
         rank = rank + digit*weight
      end do
   end function within_class_rank

   subroutine ormas_strings(space, alpha, beta, error)
      !! Every alpha and beta string of the space, in global index order
      !!
      !! When the two spins have the same electron count the two lists are
      !! identical, which the spin-transpose bookkeeping downstream depends on.
      type(ormas_space_t), intent(in) :: space
      integer(int64), allocatable, intent(out) :: alpha(:), beta(:)
      type(error_t), intent(inout) :: error

      call strings_of_one_spin(space, space%alpha_class, space%alpha_offset, alpha, error)
      if (error%has_error()) return
      call strings_of_one_spin(space, space%beta_class, space%beta_offset, beta, error)
   end subroutine ormas_strings

   pure function ormas_string_address(space, string, alpha_spin) result(index)
      !! Where a string sits in its spin's list, 1-based, or zero if it is not
      !! in the space at all
      !!
      !! Zero rather than an error: an excitation generator calls this on every
      !! string it produces, and leaving the space is the normal case ORMAS
      !! exists to express. The caller must check.
      type(ormas_space_t), intent(in) :: space
      integer(int64), intent(in) :: string
      logical, intent(in) :: alpha_spin
      integer(int64) :: index

      integer :: occupation(space%n_subspaces)
      integer :: which

      index = 0_int64
      occupation = occupation_of(space, string)

      if (alpha_spin) then
         which = class_holding(space%alpha_class, occupation)
         if (which == 0) return
         index = space%alpha_offset(which) + &
                 within_class_rank(space, space%alpha_class, which, string)
      else
         which = class_holding(space%beta_class, occupation)
         if (which == 0) return
         index = space%beta_offset(which) + &
                 within_class_rank(space, space%beta_class, which, string)
      end if
   end function ormas_string_address

   pure function one_move_apart(left, right) result(reachable)
      !! Whether two occupation vectors differ by moving one electron
      !!
      !! True for equal vectors too: an excitation inside a subspace changes no
      !! occupation at all.
      integer, intent(in) :: left(:), right(:)
      logical :: reachable

      integer :: difference(size(left))

      difference = left - right
      reachable = sum(abs(difference)) <= 2 .and. sum(difference) == 0
   end function one_move_apart

   subroutine build_ormas_closure(space, closure, error)
      !! The space widened by one excitation of either spin
      !!
      !! A pair of classes belongs if one excitation of the alpha string, or one
      !! of the beta string, could have arrived there from a pair that is in the
      !! space. Never both at once: a single excitation moves one spin.
      type(ormas_space_t), intent(in) :: space
      type(ormas_closure_t), intent(out) :: closure
      type(error_t), intent(inout) :: error

      integer :: ga, gb, source

      if (error%has_error()) return

      allocate (closure%present(space%n_beta_classes, space%n_alpha_classes))
      closure%present = .false.

      do ga = 1, space%n_alpha_classes
         do gb = 1, space%n_beta_classes
            do source = 1, space%n_alpha_classes
               if (.not. space%compatible(gb, source)) cycle
               if (one_move_apart(space%alpha_class(:, ga), space%alpha_class(:, source))) then
                  closure%present(gb, ga) = .true.
                  exit
               end if
            end do
            if (closure%present(gb, ga)) cycle
            do source = 1, space%n_beta_classes
               if (.not. space%compatible(source, ga)) cycle
               if (one_move_apart(space%beta_class(:, gb), space%beta_class(:, source))) then
                  closure%present(gb, ga) = .true.
                  exit
               end if
            end do
         end do
      end do

      call layout_addressing(closure%present, space%alpha_class_size, &
                             space%beta_class_size, space%alpha_offset, &
                             space%alpha_string_class, closure%row_length, &
                             closure%class_base, closure%block_offset, &
                             closure%alpha_base, closure%n_determinants, error)
   end subroutine build_ormas_closure

   pure function closure_address(closure, space, alpha_string, beta_string) result(address)
      !! Where a determinant sits in the widened space
      type(ormas_closure_t), intent(in) :: closure
      type(ormas_space_t), intent(in) :: space
      integer, intent(in) :: alpha_string, beta_string
      integer(int64) :: address

      integer :: ga, gb

      ga = space%alpha_string_class(alpha_string)
      gb = space%beta_string_class(beta_string)
      address = closure%alpha_base(alpha_string) + closure%block_offset(gb, ga) &
                + (int(beta_string, int64) - space%beta_offset(gb))
   end function closure_address

   pure subroutine ormas_closure_destroy(this)
      !! Give back everything the closure holds
      class(ormas_closure_t), intent(inout) :: this

      if (allocated(this%present)) deallocate (this%present)
      if (allocated(this%row_length)) deallocate (this%row_length)
      if (allocated(this%class_base)) deallocate (this%class_base)
      if (allocated(this%block_offset)) deallocate (this%block_offset)
      if (allocated(this%alpha_base)) deallocate (this%alpha_base)
      this%n_determinants = 0_int64
   end subroutine ormas_closure_destroy

   pure subroutine ormas_space_destroy(this)
      !! Give back everything the partition holds
      class(ormas_space_t), intent(inout) :: this

      if (allocated(this%first_orbital)) deallocate (this%first_orbital)
      if (allocated(this%n_orbitals)) deallocate (this%n_orbitals)
      if (allocated(this%min_electrons)) deallocate (this%min_electrons)
      if (allocated(this%max_electrons)) deallocate (this%max_electrons)
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
