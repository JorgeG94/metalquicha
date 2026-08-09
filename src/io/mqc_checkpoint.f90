!! Append-as-you-go record of computed fragment energies
module mqc_checkpoint
   !! Writes each fragment energy as it is produced, and reads them back so a
   !! rerun computes only what is missing.
   !!
   !! **The file is appended, not written at the end.** That is the whole
   !! point: a run killed at hour eleven of twelve leaves eleven hours of
   !! results behind, and the fragment table -- written once, after the
   !! expansion -- leaves nothing. Every line is flushed as it is recorded, so
   !! what survives a kill is every fragment that finished before it.
   !!
   !! **Terms are keyed by their monomers, not by their index.** A screened run
   !! and a full run number the same dimer differently, so an index-keyed
   !! checkpoint would silently pair energies with the wrong fragments the
   !! first time anyone resumed a screened list against a full one. The
   !! monomer tuple means the same thing in both.
   !!
   !! **The header carries a fingerprint and a mismatch is fatal.** Reusing an
   !! energy computed for a different geometry, basis or functional is the one
   !! failure that finishes early and reports a converged total with nothing
   !! wrong on the face of it. Refusing is the entire reason the fingerprint
   !! exists, so it is refused loudly rather than warned about.
   !!
   !! A truncated final line -- the job died mid-write -- is dropped rather
   !! than being an error. That is the expected state of a checkpoint from a
   !! killed run, which is exactly the case this is for.
   use pic_types, only: dp, int64, default_int
   use mqc_error, only: error_t, ERROR_IO, ERROR_VALIDATION
   use mqc_result_types, only: SCF_UNKNOWN
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   implicit none
   private

   public :: checkpoint_t

   character(len=*), parameter :: MAGIC = "# mqc-checkpoint v1 fingerprint="

   type :: checkpoint_t
      !! One checkpoint file, open for appending and holding what it loaded
      character(len=:), allocatable :: path
      logical :: active = .false.
      integer :: unit = -1

      integer(int64) :: n_loaded = 0
      integer :: max_level = 0
      integer, allocatable :: terms(:, :)      !! (n_loaded, max_level), sorted
      real(dp), allocatable :: energies(:)
      integer, allocatable :: scf_status(:)
   contains
      procedure :: open => checkpoint_open
      procedure :: record => checkpoint_record
      procedure :: lookup => checkpoint_lookup
      procedure :: close => checkpoint_close
   end type checkpoint_t

contains

   subroutine checkpoint_open(this, path, fingerprint, max_level, energy_only, error)
      !! Load an existing checkpoint if there is one, then open for appending
      !!
      !! Refuses rather than starting fresh when the fingerprint disagrees.
      !! Overwriting would throw away hours of a run whose file was merely
      !! misnamed; ignoring it would splice its energies into a different
      !! calculation. Neither is the caller's intent, so neither is guessed at.
      class(checkpoint_t), intent(inout) :: this
      character(len=*), intent(in) :: path
      character(len=*), intent(in) :: fingerprint
      integer, intent(in) :: max_level
      logical, intent(in) :: energy_only
         !! Whether this run computes energies alone
      type(error_t), intent(inout) :: error

      logical :: exists
      integer :: ios

      this%path = trim(path)
      this%max_level = max_level
      if (len_trim(path) == 0) return

      ! A line holds one energy and nothing else, so a reused fragment comes
      ! back without its gradient and the run dies assembling one. Refused
      ! rather than half-supported: extending the format to carry gradients
      ! and displacement sets is a real piece of work, and pretending to
      ! support them meanwhile costs a job that thought it was resumable.
      if (.not. energy_only) then
         call error%set(ERROR_VALIDATION, "checkpointing only supports energy runs; "// &
                        "this driver also needs derivatives, which the file does not "// &
                        "carry. Remove system.checkpoint, or run the energy first.")
         return
      end if

      inquire (file=this%path, exist=exists)
      if (exists) then
         call load_existing(this, fingerprint, error)
         if (error%has_error()) return
         call logger%info("Checkpoint: resuming from "//this%path//" with "// &
                          to_char(this%n_loaded)//" fragment(s) already done")
         open (newunit=this%unit, file=this%path, status="old", position="append", &
               action="write", iostat=ios)
      else
         open (newunit=this%unit, file=this%path, status="new", action="write", iostat=ios)
         if (ios == 0) then
            write (this%unit, "(a)") MAGIC//trim(fingerprint)
            flush (this%unit)
         end if
         call logger%info("Checkpoint: recording to "//this%path)
      end if

      if (ios /= 0) then
         call error%set(ERROR_IO, "could not open checkpoint "//this%path)
         return
      end if
      this%active = .true.
   end subroutine checkpoint_open

   subroutine load_existing(this, fingerprint, error)
      !! Read the header and every complete line
      class(checkpoint_t), intent(inout) :: this
      character(len=*), intent(in) :: fingerprint
      type(error_t), intent(inout) :: error

      integer :: unit, ios, level, islot, capacity
      integer(int64) :: n
      character(len=512) :: line
      character(len=:), allocatable :: stored
      integer :: row(this%max_level)
      real(dp) :: energy
      integer :: status_code

      open (newunit=unit, file=this%path, status="old", action="read", iostat=ios)
      if (ios /= 0) then
         call error%set(ERROR_IO, "could not read checkpoint "//this%path)
         return
      end if

      read (unit, "(a)", iostat=ios) line
      if (ios /= 0 .or. index(line, MAGIC) /= 1) then
         close (unit)
         call error%set(ERROR_VALIDATION, this%path//" is not a checkpoint file")
         return
      end if
      stored = trim(adjustl(line(len(MAGIC) + 1:)))
      if (stored /= trim(fingerprint)) then
         close (unit)
         call error%set(ERROR_VALIDATION, "checkpoint "//this%path//" was written by a "// &
                        "different calculation (fingerprint "//stored//", this run is "// &
                        trim(fingerprint)//"). Its energies belong to another geometry, "// &
                        "basis or method; reusing them would give a converged total for "// &
                        "neither. Point system.checkpoint at another file, or delete this one.")
         return
      end if

      capacity = 1024
      allocate (this%terms(capacity, this%max_level))
      allocate (this%energies(capacity))
      allocate (this%scf_status(capacity))
      n = 0

      do
         read (unit, "(a)", iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         ! A short read here is the last line of a job that was killed while
         ! writing it. Expected, and not an error.
         read (line, *, iostat=ios) level, (row(islot), islot=1, this%max_level), &
            energy, status_code
         if (ios /= 0) cycle
         if (level < 1 .or. level > this%max_level) cycle

         n = n + 1
         if (n > capacity) call grow(this, capacity)
         this%terms(n, :) = row
         this%energies(n) = energy
         this%scf_status(n) = status_code
      end do
      close (unit)

      this%n_loaded = n
      call sort_terms(this)
   end subroutine load_existing

   subroutine grow(this, capacity)
      !! Double the loaded arrays
      class(checkpoint_t), intent(inout) :: this
      integer, intent(inout) :: capacity

      integer, allocatable :: t(:, :), s(:)
      real(dp), allocatable :: e(:)

      allocate (t(2*capacity, this%max_level))
      allocate (e(2*capacity))
      allocate (s(2*capacity))
      t(1:capacity, :) = this%terms
      e(1:capacity) = this%energies
      s(1:capacity) = this%scf_status
      call move_alloc(t, this%terms)
      call move_alloc(e, this%energies)
      call move_alloc(s, this%scf_status)
      capacity = 2*capacity
   end subroutine grow

   subroutine checkpoint_record(this, term, energy, scf_status)
      !! Append one finished fragment, and make it survive a kill
      !!
      !! Flushed every time. A buffered checkpoint is a checkpoint that loses
      !! whatever was in the buffer, which is precisely the work a resume most
      !! wants back.
      class(checkpoint_t), intent(inout) :: this
      integer, intent(in) :: term(:)
      real(dp), intent(in) :: energy
      integer, intent(in) :: scf_status

      integer :: level, islot

      if (.not. this%active) return
      level = count(term > 0)
      write (this%unit, "(i0,*(1x,i0))", advance="no") level, &
         (term(islot), islot=1, this%max_level)
      write (this%unit, "(1x,es24.16,1x,i0)") energy, scf_status
      flush (this%unit)
   end subroutine checkpoint_record

   subroutine checkpoint_lookup(this, term, found, energy, scf_status)
      !! Is this term already done, and with what energy
      !!
      !! Bisection over the sorted load rather than a scan: a resume asks this
      !! once per term, and at the fragment counts this exists for, linear
      !! search would cost more than recomputing everything.
      class(checkpoint_t), intent(in) :: this
      integer, intent(in) :: term(:)
      logical, intent(out) :: found
      real(dp), intent(out) :: energy
      integer, intent(out) :: scf_status

      integer(int64) :: lo, hi, mid
      integer :: order

      found = .false.
      energy = 0.0_dp
      scf_status = SCF_UNKNOWN
      if (this%n_loaded <= 0) return

      lo = 1
      hi = this%n_loaded
      do while (lo <= hi)
         mid = (lo + hi)/2
         order = compare(this%terms(mid, :), term, this%max_level)
         if (order == 0) then
            found = .true.
            energy = this%energies(mid)
            scf_status = this%scf_status(mid)
            return
         end if
         if (order < 0) then
            lo = mid + 1
         else
            hi = mid - 1
         end if
      end do
   end subroutine checkpoint_lookup

   pure function compare(a, b, n) result(order)
      !! Lexicographic order on two zero-padded monomer rows
      integer, intent(in) :: a(:), b(:)
      integer, intent(in) :: n
      integer :: order

      integer :: i

      order = 0
      do i = 1, n
         if (a(i) < b(i)) then
            order = -1
            return
         end if
         if (a(i) > b(i)) then
            order = 1
            return
         end if
      end do
   end function compare

   subroutine sort_terms(this)
      !! Insertion sort by monomer tuple
      !!
      !! Quadratic, and deliberately so for now: the load is bounded by what a
      !! previous run finished, and anything large enough for this to hurt is
      !! also large enough that the sort wants to happen while reading rather
      !! than after. Marked so the next person knows it was a choice.
      class(checkpoint_t), intent(inout) :: this

      integer(int64) :: i, j
      integer :: key_term(this%max_level)
      integer :: key_status
      real(dp) :: key_energy

      do i = 2, this%n_loaded
         key_term = this%terms(i, :)
         key_energy = this%energies(i)
         key_status = this%scf_status(i)
         j = i - 1
         do while (j >= 1)
            if (compare(this%terms(j, :), key_term, this%max_level) <= 0) exit
            this%terms(j + 1, :) = this%terms(j, :)
            this%energies(j + 1) = this%energies(j)
            this%scf_status(j + 1) = this%scf_status(j)
            j = j - 1
         end do
         this%terms(j + 1, :) = key_term
         this%energies(j + 1) = key_energy
         this%scf_status(j + 1) = key_status
      end do
   end subroutine sort_terms

   subroutine checkpoint_close(this)
      !! Close the file. Loaded entries stay readable.
      class(checkpoint_t), intent(inout) :: this

      if (this%active .and. this%unit >= 0) close (this%unit)
      this%active = .false.
      this%unit = -1
   end subroutine checkpoint_close

end module mqc_checkpoint
