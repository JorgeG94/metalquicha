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
   use mqc_io_helpers, only: ends_with
   use mqc_hdf5_checkpoint, only: hdf5_checkpoint_t, hdf5_checkpoint_available
   implicit none
   private

   public :: checkpoint_t

   character(len=*), parameter :: MAGIC = "# mqc-checkpoint v1 fingerprint="

   type :: checkpoint_t
      !! One checkpoint file, open for appending and holding what it loaded
      character(len=:), allocatable :: path
      logical :: active = .false.
      integer :: unit = 0
         !! Meaningful only while `active`; `newunit=` values are negative,
         !! so the sign of this says nothing about whether it is open.

      integer(int64) :: n_loaded = 0
      integer :: max_level = 0
      integer, allocatable :: terms(:, :)      !! (n_loaded, max_level), sorted
      real(dp), allocatable :: energies(:)
      integer, allocatable :: scf_status(:)
      real(dp), allocatable :: homo(:), lumo(:)
      logical, allocatable :: has_orbitals(:)

      logical :: use_hdf5 = .false.
         !! Which backend is in use. Text holds one energy per record and
         !! needs no dependency, which is right for the screening passes that
         !! are the common case; HDF5 holds derivatives as well and is chosen
         !! when the driver needs them, or when the path asks for it by name.
      type(hdf5_checkpoint_t) :: h5
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
      ! The backend is chosen, not configured. A run needing derivatives has
      ! to have HDF5 because the text format holds one energy per record; an
      ! energy-only run gets text because it is crash-proof line by line and
      ! costs no dependency. A caller can still ask for HDF5 by naming the
      ! file .h5, which is what a screening pass wants when it is large.
      this%use_hdf5 = (.not. energy_only) .or. ends_with(this%path, ".h5") &
                      .or. ends_with(this%path, ".hdf5")

      if (this%use_hdf5 .and. .not. hdf5_checkpoint_available()) then
         call error%set(ERROR_VALIDATION, "this run needs an HDF5 checkpoint -- a "// &
                        "derivative run, or a .h5 path -- and this build has none. "// &
                        "Configure with -DMQC_ENABLE_HDF5=ON, or run the energy "// &
                        "first and checkpoint that.")
         return
      end if

      if (this%use_hdf5) then
         call this%h5%open(this%path, fingerprint, max_level, error)
         if (error%has_error()) return
         this%active = .true.
         call logger%info("Checkpoint: HDF5 at "//this%path//" ("// &
                          to_char(this%h5%n_loaded)//" fragment(s) already done)")
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
      real(dp) :: homo_in, lumo_in
      integer :: orbitals_in

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
      allocate (this%homo(capacity), this%lumo(capacity), this%has_orbitals(capacity))
      n = 0

      do
         read (unit, "(a)", iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         ! A short read here is the last line of a job that was killed while
         ! writing it. Expected, and not an error.
         ! Orbitals are read separately so a checkpoint written before this
         ! column existed still loads -- it simply has no frontier pair.
         orbitals_in = 0
         homo_in = 0.0_dp
         lumo_in = 0.0_dp
         read (line, *, iostat=ios) level, (row(islot), islot=1, this%max_level), &
            energy, status_code, orbitals_in, homo_in, lumo_in
         if (ios /= 0) then
            orbitals_in = 0
            read (line, *, iostat=ios) level, (row(islot), islot=1, this%max_level), &
               energy, status_code
         end if
         if (ios /= 0) cycle
         if (level < 1 .or. level > this%max_level) cycle

         n = n + 1
         if (n > capacity) call grow(this, capacity)
         this%terms(n, :) = row
         this%energies(n) = energy
         this%scf_status(n) = status_code
         this%homo(n) = homo_in
         this%lumo(n) = lumo_in
         this%has_orbitals(n) = (orbitals_in /= 0)
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
      block
         real(dp), allocatable :: h(:), l(:)
         logical, allocatable :: o(:)
         allocate (h(2*capacity), l(2*capacity), o(2*capacity))
         h(1:capacity) = this%homo
         l(1:capacity) = this%lumo
         o(1:capacity) = this%has_orbitals
         call move_alloc(h, this%homo)
         call move_alloc(l, this%lumo)
         call move_alloc(o, this%has_orbitals)
      end block
      capacity = 2*capacity
   end subroutine grow

   subroutine checkpoint_record(this, term, energy, scf_status, n_atoms, gradient, hessian, &
                                homo, lumo, has_orbitals)
      !! Append one finished fragment, and make it survive a kill
      !!
      !! Flushed every time. A buffered checkpoint is a checkpoint that loses
      !! whatever was in the buffer, which is precisely the work a resume most
      !! wants back.
      class(checkpoint_t), intent(inout) :: this
      integer, intent(in) :: term(:)
      real(dp), intent(in) :: energy
      integer, intent(in) :: scf_status
      integer, intent(in), optional :: n_atoms
         !! Atoms in the fragment, caps included. Needed to give a stored
         !! gradient its shape back; ignored by the text backend.
      real(dp), intent(in), optional :: gradient(:, :)
      real(dp), intent(in), optional :: hessian(:, :)
      real(dp), intent(in), optional :: homo, lumo
      logical, intent(in), optional :: has_orbitals
         !! Whether the method reported a frontier pair. Passed as a flag
         !! rather than by omitting the values, so a caller forwards three
         !! arguments unconditionally instead of branching -- and a resumed
         !! fragment does not come back claiming a gap of zero.

      integer :: level, islot, natoms_local
      logical :: orbitals_out

      if (.not. this%active) return

      if (this%use_hdf5) then
         natoms_local = 0
         if (present(n_atoms)) natoms_local = n_atoms
         if (present(gradient) .and. present(hessian)) then
            call this%h5%record(term, energy, scf_status, natoms_local, gradient, hessian)
         else if (present(gradient)) then
            call this%h5%record(term, energy, scf_status, natoms_local, gradient=gradient)
         else if (present(hessian)) then
            call this%h5%record(term, energy, scf_status, natoms_local, hessian=hessian)
         else
            call this%h5%record(term, energy, scf_status, natoms_local)
         end if
         return
      end if
      level = count(term > 0)
      write (this%unit, "(i0,*(1x,i0))", advance="no") level, &
         (term(islot), islot=1, this%max_level)
      write (this%unit, "(1x,es24.16,1x,i0)", advance="no") energy, scf_status
      orbitals_out = .false.
      if (present(has_orbitals)) orbitals_out = has_orbitals
      if (orbitals_out .and. present(homo) .and. present(lumo)) then
         write (this%unit, "(1x,i0,2(1x,es24.16))") 1, homo, lumo
      else
         write (this%unit, "(1x,i0,2(1x,es24.16))") 0, 0.0_dp, 0.0_dp
      end if
      flush (this%unit)
   end subroutine checkpoint_record

   subroutine checkpoint_lookup(this, term, found, energy, scf_status, n_atoms, gradient, &
                                hessian, homo, lumo, has_orbitals)
      !! Is this term already done, and with what energy
      !!
      !! Bisection over the sorted load rather than a scan: a resume asks this
      !! once per term, and at the fragment counts this exists for, linear
      !! search would cost more than recomputing everything.
      class(checkpoint_t), intent(inout) :: this
      integer, intent(in) :: term(:)
      logical, intent(out) :: found
      real(dp), intent(out) :: energy
      integer, intent(out) :: scf_status
      integer, intent(out), optional :: n_atoms
      real(dp), allocatable, intent(out), optional :: gradient(:, :)
      real(dp), allocatable, intent(out), optional :: hessian(:, :)
      real(dp), intent(out), optional :: homo, lumo
      logical, intent(out), optional :: has_orbitals

      integer(int64) :: lo, hi, mid
      integer :: order, natoms_local

      found = .false.
      energy = 0.0_dp
      scf_status = SCF_UNKNOWN
      if (present(n_atoms)) n_atoms = 0
      if (present(homo)) homo = 0.0_dp
      if (present(lumo)) lumo = 0.0_dp
      if (present(has_orbitals)) has_orbitals = .false.

      if (this%use_hdf5) then
         if (present(gradient) .and. present(hessian)) then
            call this%h5%lookup(term, found, energy, scf_status, natoms_local, &
                                gradient=gradient, hessian=hessian)
         else if (present(gradient)) then
            call this%h5%lookup(term, found, energy, scf_status, natoms_local, gradient=gradient)
         else if (present(hessian)) then
            call this%h5%lookup(term, found, energy, scf_status, natoms_local, hessian=hessian)
         else
            call this%h5%lookup(term, found, energy, scf_status, natoms_local)
         end if
         if (present(n_atoms)) n_atoms = natoms_local
         return
      end if

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
            if (present(homo)) homo = this%homo(mid)
            if (present(lumo)) lumo = this%lumo(mid)
            if (present(has_orbitals)) has_orbitals = this%has_orbitals(mid)
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
      real(dp) :: key_homo, key_lumo
      logical :: key_orbitals

      do i = 2, this%n_loaded
         key_term = this%terms(i, :)
         key_energy = this%energies(i)
         key_status = this%scf_status(i)
         key_homo = this%homo(i)
         key_lumo = this%lumo(i)
         key_orbitals = this%has_orbitals(i)
         j = i - 1
         do while (j >= 1)
            if (compare(this%terms(j, :), key_term, this%max_level) <= 0) exit
            this%terms(j + 1, :) = this%terms(j, :)
            this%energies(j + 1) = this%energies(j)
            this%scf_status(j + 1) = this%scf_status(j)
            this%homo(j + 1) = this%homo(j)
            this%lumo(j + 1) = this%lumo(j)
            this%has_orbitals(j + 1) = this%has_orbitals(j)
            j = j - 1
         end do
         this%terms(j + 1, :) = key_term
         this%energies(j + 1) = key_energy
         this%scf_status(j + 1) = key_status
         this%homo(j + 1) = key_homo
         this%lumo(j + 1) = key_lumo
         this%has_orbitals(j + 1) = key_orbitals
      end do
   end subroutine sort_terms

   subroutine checkpoint_close(this)
      !! Close the file. Loaded entries stay readable.
      class(checkpoint_t), intent(inout) :: this

      if (this%use_hdf5) then
         call this%h5%close()
      else if (this%active) then
         ! Not `this%unit >= 0`: `newunit=` hands back *negative* unit
         ! numbers, so that test was never true and the file was never
         ! closed. A program that runs once and exits does not notice --
         ! the process closes it -- but a session runs many calculations
         ! in one process, and the second one asking for the same
         ! checkpoint got "already opened on another unit" instead of a
         ! resume. `active` is the thing that actually says a unit is open.
         close (this%unit)
      end if
      this%active = .false.
      this%unit = 0
   end subroutine checkpoint_close

end module mqc_checkpoint
