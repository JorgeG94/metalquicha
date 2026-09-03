!! Fragment term list as a standalone, inspectable object
module mqc_fraglist
   !! A list of MBE terms that can be built, handed out, filtered and handed
   !! back, rather than generated and consumed inside one call.
   !!
   !! A screen the MBE engine applies on the way past can only ever be a
   !! distance screen. Keeping terms whose `delta_energy` from a previous run
   !! cleared a threshold, resuming a run that died part way, or hand-picking a
   !! set all need the list to exist as something a caller can hold.
   !!
   !! Layout. Terms are **monomer index sets**, not atom sets, and the store is
   !! rectangular `(n_terms, max_level)` with unused slots zero -- one term per
   !! ROW, which is the order `generate_fragment_list` writes and every MBE
   !! consumer reads (`count(polymers(i, :) > 0)` is the order of term i). It
   !! is also the shape the per-fragment CSV writes, `m1..m<max_level>`
   !! zero-filled, so a term here and a row there line up.
   !!
   !! Indices are 1-based monomer numbers, matching `create_monomer_list` and
   !! the CSV. That differs from the 0-based *atom* indices the input schema
   !! uses, which is a real trap: these are different things counted
   !! differently, and the CSV is the one to check against.
   use pic_types, only: dp, default_int, int64
   use mqc_combinatorics, only: get_nfrags, create_monomer_list, generate_fragment_list, &
                                calculate_fragment_distances
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: fraglist_t

   integer, parameter :: MAX_TERM_LEVEL = 10
      !! Deepest expansion subset closure handles, matching MAX_MBE_LEVEL. The
      !! closure works on stack arrays of this width to stay out of the heap in
      !! a loop that runs once per term.

   type :: fraglist_t
      !! A term list, and whatever has been computed about its terms
      integer(default_int) :: max_level = 0
         !! Columns in `terms`; a term of lower order zero-fills the rest
      integer(int64) :: n_terms = 0
      integer(default_int), allocatable :: terms(:, :)
         !! (n_terms, max_level), 1-based monomer indices, zero-padded
      real(dp), allocatable :: distances(:)
         !! Minimum inter-monomer distance per term, Angstrom. Allocated only
         !! once `measure` has run -- a list that was set by a caller has no
         !! distances until it is asked for them.
      logical :: has_distances = .false.
   contains
      procedure :: create => fraglist_create
      procedure :: replace => fraglist_replace
      procedure :: measure => fraglist_measure
      procedure :: keep => fraglist_keep
      procedure :: close_subsets => fraglist_close_subsets
      procedure :: destroy => fraglist_destroy
      procedure :: level_of => fraglist_level_of
   end type fraglist_t

contains

   subroutine fraglist_create(this, n_monomers, max_level, error)
      !! Every combination of monomers from 2 up to max_level
      !!
      !! Monomers themselves are not terms here. The MBE evaluates them, but
      !! they are implied by the monomer count and `generate_fragment_list`
      !! starts at pairs; a caller that wants them in the list can add them.
      class(fraglist_t), intent(inout) :: this
      integer(default_int), intent(in) :: n_monomers
      integer(default_int), intent(in) :: max_level
      type(error_t), intent(inout) :: error

      integer(default_int), allocatable :: monomers(:)
      integer(int64) :: capacity, count

      call this%destroy()

      if (n_monomers <= 0) then
         call error%set(ERROR_VALIDATION, "fraglist: a system needs at least one monomer")
         return
      end if
      if (max_level < 2) then
         call error%set(ERROR_VALIDATION, &
                        "fraglist: max_level must be at least 2 to have any n-mer terms")
         return
      end if
      if (max_level > n_monomers) then
         call error%set(ERROR_VALIDATION, &
                        "fraglist: max_level exceeds the number of monomers")
         return
      end if

      ! get_nfrags counts monomers too, so it is an upper bound rather than the
      ! exact size -- generous by n_monomers entries, which is nothing beside
      ! the n-mer count and saves counting the combinations twice.
      capacity = get_nfrags(n_monomers, max_level)

      allocate (monomers(n_monomers))
      call create_monomer_list(monomers)

      allocate (this%terms(capacity, max_level))
      this%terms = 0
      count = 0_int64
      call generate_fragment_list(monomers, max_level, this%terms, count)
      deallocate (monomers)

      this%max_level = max_level
      this%n_terms = count
   end subroutine fraglist_create

   subroutine fraglist_replace(this, terms, n_terms, max_level, error)
      !! Take a caller's term list wholesale
      !!
      !! Whatever was known about the old list -- distances above all -- does
      !! not describe the new one, so it is dropped rather than left to be read
      !! as though it did.
      class(fraglist_t), intent(inout) :: this
      integer(default_int), intent(in) :: terms(:, :)  !! (n_terms, max_level)
      integer(int64), intent(in) :: n_terms
      integer(default_int), intent(in) :: max_level
      type(error_t), intent(inout) :: error

      if (n_terms < 0) then
         call error%set(ERROR_VALIDATION, "fraglist: term count cannot be negative")
         return
      end if
      if (max_level < 1) then
         call error%set(ERROR_VALIDATION, "fraglist: max_level must be at least 1")
         return
      end if
      if (size(terms, 1) < n_terms .or. size(terms, 2) < max_level) then
         call error%set(ERROR_VALIDATION, &
                        "fraglist: the supplied array is smaller than the term count claims")
         return
      end if

      call this%destroy()
      this%max_level = max_level
      this%n_terms = n_terms
      allocate (this%terms(max(n_terms, 1_int64), max_level))
      this%terms = 0
      if (n_terms > 0) this%terms(1:n_terms, 1:max_level) = terms(1:n_terms, 1:max_level)
   end subroutine fraglist_replace

   subroutine fraglist_measure(this, sys_geom, error)
      !! Compute the minimum inter-monomer distance of every term
      !!
      !! Separate from `create` because a list that arrived from a caller needs
      !! it too, and because it is `O(n_terms * atoms^2)` and a caller screening
      !! on energy never needs it at all.
      class(fraglist_t), intent(inout) :: this
      type(system_geometry_t), intent(in) :: sys_geom
      type(error_t), intent(inout) :: error

      if (this%n_terms <= 0) then
         call error%set(ERROR_VALIDATION, "fraglist: nothing to measure, the list is empty")
         return
      end if

      if (allocated(this%distances)) deallocate (this%distances)
      allocate (this%distances(this%n_terms))
      call calculate_fragment_distances(this%terms, this%n_terms, sys_geom, this%distances)
      this%has_distances = .true.
   end subroutine fraglist_measure

   subroutine fraglist_keep(this, mask, error)
      !! Keep the terms the mask selects, in their existing order
      !!
      !! The one operation every screen reduces to: distance screening, an
      !! energy threshold read out of a previous run's CSV, and resuming from
      !! the terms a dead run never reached are this call with a different mask.
      class(fraglist_t), intent(inout) :: this
      logical, intent(in) :: mask(:)  !! One entry per term
      type(error_t), intent(inout) :: error

      integer(default_int), allocatable :: kept(:, :)
      real(dp), allocatable :: kept_distances(:)
      integer(int64) :: iterm, n_kept

      if (size(mask, kind=int64) /= this%n_terms) then
         call error%set(ERROR_VALIDATION, &
                        "fraglist: the mask does not have one entry per term")
         return
      end if

      n_kept = count(mask, kind=int64)
      allocate (kept(max(n_kept, 1_int64), this%max_level))
      kept = 0
      if (this%has_distances) allocate (kept_distances(max(n_kept, 1_int64)))

      n_kept = 0
      do iterm = 1, this%n_terms
         if (.not. mask(iterm)) cycle
         n_kept = n_kept + 1
         kept(n_kept, :) = this%terms(iterm, :)
         if (this%has_distances) kept_distances(n_kept) = this%distances(iterm)
      end do

      call move_alloc(kept, this%terms)
      if (this%has_distances) call move_alloc(kept_distances, this%distances)
      this%n_terms = n_kept
   end subroutine fraglist_keep

   subroutine fraglist_close_subsets(this, error)
      !! Add back every subset of every surviving term
      !!
      !! **A screened list is not usable until this has run.** An n-body term's
      !! contribution is its energy less the delta of every proper subset, so
      !! keeping a trimer whose dimers fell below the threshold does not
      !! approximate anything -- the assembly looks the subsets up and fails.
      !! Screening therefore selects which *high-order* terms to keep, and this
      !! restores what computing them requires.
      !!
      !! That is also why an energy screen saves less than it first appears: a
      !! kept trimer drags in three dimers and three monomers, most of which
      !! the screen had already rejected.
      ! TODO(mqc): quadratic, and on the path every screened run has to take.
      ! `add_unique` scans everything already stored for each of the 2^level
      ! subsets of every term, so the cost grows with the square of the closed
      ! list. The store is also preallocated at `n_terms * 2**max_level` rows,
      ! eight times the term count at level 3, before any of it is filled.
      class(fraglist_t), intent(inout) :: this
      type(error_t), intent(inout) :: error

      integer(default_int), allocatable :: closed(:, :)
      integer(default_int) :: subset(MAX_TERM_LEVEL), term(MAX_TERM_LEVEL)
      integer(default_int) :: pick(MAX_TERM_LEVEL)
      integer(int64) :: iterm, n_closed, capacity
      integer(default_int) :: level, size_wanted, i
      logical :: has_next

      if (this%n_terms <= 0) return
      if (this%max_level > MAX_TERM_LEVEL) then
         call error%set(ERROR_VALIDATION, &
                        "fraglist: max_level exceeds what subset closure supports")
         return
      end if

      ! Every term of order n contributes at most 2^n - 1 subsets including
      ! itself; duplicates are dropped below, so this is a bound, not a count.
      capacity = this%n_terms*int(2**this%max_level, int64)
      allocate (closed(capacity, this%max_level))
      closed = 0
      n_closed = 0

      do iterm = 1, this%n_terms
         term = 0
         level = this%level_of(iterm)
         term(1:level) = this%terms(iterm, 1:level)
         do size_wanted = 1, level
            call first_combination(subset, size_wanted)
            do
               ! Gathered into a local rather than passed as
               ! `[(term(subset(i)), i=1, size_wanted)]`. An array constructor
               ! as an actual argument builds a temporary, and nvfortran sizes
               ! that temporary wrongly under `-Mvect=simd`, which `-fast`
               ! implies: the run aborts asking for a garbage allocation whose
               ! size varies run to run. `-O2` alone and `-fast -Mnovect` are
               ! both fine.
               do i = 1, size_wanted
                  pick(i) = term(subset(i))
               end do
               call add_unique(closed, n_closed, this%max_level, &
                               pick(1:size_wanted), size_wanted)
               call next_index_set(subset, size_wanted, level, has_next)
               if (.not. has_next) exit
            end do
         end do
      end do

      if (allocated(this%distances)) deallocate (this%distances)
      this%has_distances = .false.
      call move_alloc(closed, this%terms)
      this%n_terms = n_closed
   end subroutine fraglist_close_subsets

   subroutine add_unique(store, n_stored, max_level, candidate, level)
      !! Append a term unless an identical one is already present
      integer(default_int), intent(inout) :: store(:, :)
      integer(int64), intent(inout) :: n_stored
      integer(default_int), intent(in) :: max_level
      integer(default_int), intent(in) :: candidate(:)
      integer(default_int), intent(in) :: level

      integer(int64) :: istored

      do istored = 1, n_stored
         if (count(store(istored, :) > 0) /= level) cycle
         if (all(store(istored, 1:level) == candidate(1:level))) return
      end do

      n_stored = n_stored + 1
      store(n_stored, :) = 0
      store(n_stored, 1:level) = candidate(1:level)
   end subroutine add_unique

   pure subroutine first_combination(indices, k)
      !! The lexicographically first choice of k positions
      integer(default_int), intent(inout) :: indices(:)
      integer(default_int), intent(in) :: k
      integer(default_int) :: i

      do i = 1, k
         indices(i) = i
      end do
   end subroutine first_combination

   pure subroutine next_index_set(indices, k, n, has_next)
      !! Advance k chosen positions out of n, lexicographically
      integer(default_int), intent(inout) :: indices(:)
      integer(default_int), intent(in) :: k, n
      logical, intent(out) :: has_next
      integer(default_int) :: i, j

      has_next = .false.
      do i = k, 1, -1
         if (indices(i) < n - (k - i)) then
            indices(i) = indices(i) + 1
            do j = i + 1, k
               indices(j) = indices(j - 1) + 1
            end do
            has_next = .true.
            return
         end if
      end do
   end subroutine next_index_set

   pure function fraglist_level_of(this, iterm) result(level)
      !! How many monomers a term actually has
      !!
      !! The store is padded to `max_level`, so the order of a term is the
      !! count of its non-zero entries rather than the number of rows.
      class(fraglist_t), intent(in) :: this
      integer(int64), intent(in) :: iterm
      integer(default_int) :: level

      level = 0
      if (iterm < 1 .or. iterm > this%n_terms) return
      level = count(this%terms(iterm, :) > 0)
   end function fraglist_level_of

   subroutine fraglist_destroy(this)
      !! Release the term store and anything measured about it
      class(fraglist_t), intent(inout) :: this

      if (allocated(this%terms)) deallocate (this%terms)
      if (allocated(this%distances)) deallocate (this%distances)
      this%max_level = 0
      this%n_terms = 0
      this%has_distances = .false.
   end subroutine fraglist_destroy

end module mqc_fraglist
