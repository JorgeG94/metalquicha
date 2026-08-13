!! Read a `.efp` fragment potential back in
module mqc_efp_read
   !! The other half of `mqc_efp_potential`: that module computes a potential and
   !! writes it, this one reads one and hands back the parameters an interaction
   !! energy needs.
   !!
   !! **Only the sections an energy reads are parsed.** A potential carries
   !! seventeen sections and most of them are bulky -- the dynamic polarizabilities
   !! alone are twelve frequencies by however many localized orbitals, and the
   !! projection data carries a whole basis set. Electrostatics needs the expansion
   !! points, the multipoles and the two screening fits, and that is what this reads.
   !! The rest is skipped by name rather than by position, so adding a section later
   !! is a matter of naming it.
   !!
   !! **The format is GAMESS's, so it is read the way GAMESS writes it** rather than
   !! the way a reader would prefer:
   !!
   !!   * a section is a header line, records, then `STOP`
   !!   * a record is a label followed by numbers, and a record too long for a line
   !!     is continued with a trailing `>`
   !!   * a header is distinguished from a record by starting in column two and
   !!     carrying no numbers of its own -- except that some do (`MULTIPLICITY 1`,
   !!     `PROJECTION WAVEFUNCTION 4 19`), which is why the known names are matched
   !!     rather than the shape guessed at
   !!
   !! Units are as the file carries them: Bohr for the points, atomic units for the
   !! multipoles. Nothing is converted on the way in.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION, ERROR_IO
   implicit none
   private

   public :: efp_fragment_t
   public :: read_efp_potential

   !> Multipole components at one expansion point, in the file's own order.
   integer, parameter :: N_DIPOLE = 3
   integer, parameter :: N_QUADRUPOLE = 6
   integer, parameter :: N_OCTUPOLE = 10

   !> A dynamic polarizability record: a centroid then a 3x3 tensor.
   integer, parameter :: N_DYNAMIC_RECORD = 12

   !> Row and column of each of GAMESS's nine polarizability slots. The diagonal
   !> comes first and the off-diagonal triples are the transpose of what its labels
   !> suggest -- the same map `mqc_efp_potential` writes with, established there
   !> against GAMESS's own output. Reading these nine as a row-major 3x3 gives a
   !> tensor whose trace is negative, which is how the mistake announces itself.
   integer, parameter :: N_POL_SLOTS = 9
   integer, parameter :: POL_ROW(N_POL_SLOTS) = [1, 2, 3, 2, 3, 3, 1, 1, 2]
   integer, parameter :: POL_COL(N_POL_SLOTS) = [1, 2, 3, 1, 1, 2, 2, 3, 3]

   !> Longest line a potential is expected to carry, matching the writer's own.
   integer, parameter :: MAX_LINE = 160

   type :: efp_fragment_t
      !! One fragment's parameters, as read
      character(len=:), allocatable :: name        !! From `$FRAGNAME`, without the `$`
      integer :: n_points = 0                      !! Atoms, then bond midpoints
      integer :: n_atoms = 0                       !! Points carrying a nuclear charge
      integer :: multiplicity = 1
      character(len=8), allocatable :: labels(:)
      real(dp), allocatable :: points(:, :)        !! (3, n_points), Bohr
      real(dp), allocatable :: mass(:)             !! amu, zero at a midpoint
      real(dp), allocatable :: charge(:)           !! Z, zero at a midpoint
      real(dp), allocatable :: q_elec(:)           !! Electronic monopole
      real(dp), allocatable :: q_nuc(:)            !! Nuclear monopole
      real(dp), allocatable :: dipole(:, :)        !! (3, n_points)
      real(dp), allocatable :: quadrupole(:, :)    !! (6, n_points)
      real(dp), allocatable :: octopole(:, :)      !! (10, n_points)
      real(dp), allocatable :: screen(:)           !! Gaussian damping exponent per point
      real(dp), allocatable :: screen2(:)          !! Exponential damping exponent per point
      logical :: has_screen = .false.
      logical :: has_screen2 = .false.
      !> Dynamic polarizabilities, for dispersion. `(3, 3, n_lmo, n_freq)`.
      integer :: n_lmo = 0
      integer :: n_freq = 0
      real(dp), allocatable :: dyn_pol(:, :, :, :)
      real(dp), allocatable :: centroids(:, :)     !! (3, n_lmo), Bohr
      real(dp), allocatable :: frequencies(:)      !! Imaginary, a.u.
      logical :: has_dynamic = .false.
      !> Static polarizabilities at the same centroids, for polarization.
      integer :: n_pol = 0
      real(dp), allocatable :: static_pol(:, :, :)  !! (3, 3, n_pol)
      real(dp), allocatable :: pol_points(:, :)     !! (3, n_pol), Bohr
      logical :: has_static_pol = .false.
   contains
      procedure :: destroy => fragment_destroy
      procedure :: net_charge => fragment_net_charge
   end type efp_fragment_t

contains

   subroutine fragment_destroy(self)
      class(efp_fragment_t), intent(inout) :: self

      if (allocated(self%name)) deallocate (self%name)
      if (allocated(self%labels)) deallocate (self%labels)
      if (allocated(self%points)) deallocate (self%points)
      if (allocated(self%mass)) deallocate (self%mass)
      if (allocated(self%charge)) deallocate (self%charge)
      if (allocated(self%q_elec)) deallocate (self%q_elec)
      if (allocated(self%q_nuc)) deallocate (self%q_nuc)
      if (allocated(self%dipole)) deallocate (self%dipole)
      if (allocated(self%quadrupole)) deallocate (self%quadrupole)
      if (allocated(self%octopole)) deallocate (self%octopole)
      if (allocated(self%screen)) deallocate (self%screen)
      if (allocated(self%screen2)) deallocate (self%screen2)
      if (allocated(self%dyn_pol)) deallocate (self%dyn_pol)
      if (allocated(self%centroids)) deallocate (self%centroids)
      if (allocated(self%frequencies)) deallocate (self%frequencies)
      self%n_lmo = 0
      self%n_freq = 0
      self%has_dynamic = .false.
      if (allocated(self%static_pol)) deallocate (self%static_pol)
      if (allocated(self%pol_points)) deallocate (self%pol_points)
      self%n_pol = 0
      self%has_static_pol = .false.
      self%n_points = 0
      self%n_atoms = 0
      self%has_screen = .false.
      self%has_screen2 = .false.
   end subroutine fragment_destroy

   pure function fragment_net_charge(self) result(q)
      !! The fragment's net charge, nuclear plus electronic
      !!
      !! Worth having as its own function because it is the cheapest check that a
      !! potential was read correctly: the monopoles must sum to the charge the
      !! fragment was computed for, and a misparsed record shows up here long before
      !! it shows up in an energy.
      class(efp_fragment_t), intent(in) :: self
      real(dp) :: q

      q = 0.0_dp
      if (allocated(self%q_elec)) q = q + sum(self%q_elec)
      if (allocated(self%q_nuc)) q = q + sum(self%q_nuc)
   end function fragment_net_charge

   subroutine read_efp_potential(path, frag, error)
      !! Read the electrostatic parameters out of a `.efp` file
      character(len=*), intent(in) :: path
      type(efp_fragment_t), intent(out) :: frag
      type(error_t), intent(inout) :: error

      character(len=MAX_LINE), allocatable :: lines(:)
      character(len=MAX_LINE), allocatable :: labels(:)
      real(dp), allocatable :: values(:, :)
      integer :: n_lines, i, n_rec
      character(len=:), allocatable :: section

      call slurp(path, lines, n_lines, error)
      if (error%has_error()) return

      call fragment_name(lines, n_lines, frag%name)

      ! COORDINATES first: it fixes how many points every other section must carry,
      ! so a section that disagrees is caught rather than silently truncated.
      call section_records(lines, n_lines, "COORDINATES", 5, labels, values, n_rec, error)
      if (error%has_error()) return
      if (n_rec == 0) then
         call error%set(ERROR_VALIDATION, "efp: '"//trim(path)// &
                        "' carries no COORDINATES section")
         return
      end if
      frag%n_points = n_rec
      allocate (frag%labels(n_rec), frag%points(3, n_rec), frag%mass(n_rec), &
                frag%charge(n_rec))
      do i = 1, n_rec
         frag%labels(i) = trim(labels(i))
         frag%points(:, i) = values(1:3, i)
         frag%mass(i) = values(4, i)
         frag%charge(i) = values(5, i)
      end do
      frag%n_atoms = count(frag%charge > 0.0_dp)

      call one_section(lines, n_lines, "MONOPOLES", 2, frag%n_points, values, error)
      if (error%has_error()) return
      allocate (frag%q_elec(frag%n_points), frag%q_nuc(frag%n_points))
      frag%q_elec = values(1, :)
      frag%q_nuc = values(2, :)

      call one_section(lines, n_lines, "DIPOLES", N_DIPOLE, frag%n_points, values, error)
      if (error%has_error()) return
      allocate (frag%dipole(N_DIPOLE, frag%n_points))
      frag%dipole = values

      call one_section(lines, n_lines, "QUADRUPOLES", N_QUADRUPOLE, frag%n_points, &
                       values, error)
      if (error%has_error()) return
      allocate (frag%quadrupole(N_QUADRUPOLE, frag%n_points))
      frag%quadrupole = values

      call one_section(lines, n_lines, "OCTUPOLES", N_OCTUPOLE, frag%n_points, &
                       values, error)
      if (error%has_error()) return
      allocate (frag%octopole(N_OCTUPOLE, frag%n_points))
      frag%octopole = values

      ! The two screening sections are optional: a potential written without a
      ! damping fit is still a usable potential, it just has no penetration term.
      ! Each record is `1.0 alpha`, the leading one being a weight GAMESS always
      ! writes as unity.
      call section_records(lines, n_lines, "SCREEN2", 2, labels, values, n_rec, error)
      if (error%has_error()) return
      if (n_rec > 0) then
         call expect_points(n_rec, frag%n_points, "SCREEN2", error)
         if (error%has_error()) return
         allocate (frag%screen2(frag%n_points))
         frag%screen2 = values(2, :)
         frag%has_screen2 = .true.
      end if

      call section_records(lines, n_lines, "SCREEN", 2, labels, values, n_rec, error)
      if (error%has_error()) return
      if (n_rec > 0) then
         call expect_points(n_rec, frag%n_points, "SCREEN", error)
         if (error%has_error()) return
         allocate (frag%screen(frag%n_points))
         frag%screen = values(2, :)
         frag%has_screen = .true.
      end if

      call read_static_pol(lines, n_lines, frag, error)
      if (error%has_error()) return

      call read_dynamic(lines, n_lines, frag, error)
      if (error%has_error()) return

      call multiplicity(lines, n_lines, frag%multiplicity)
      deallocate (lines)
      if (allocated(labels)) deallocate (labels)
      if (allocated(values)) deallocate (values)
   end subroutine read_efp_potential

   subroutine read_static_pol(lines, n_lines, frag, error)
      !! `POLARIZABLE POINTS`: the static polarizability at each orbital centroid
      !!
      !! The same record shape as the dynamic section -- a label, a centroid, nine
      !! tensor components in GAMESS's slot order -- with one block rather than
      !! twelve, and **one** label token rather than two: `CT1` here against `CT  1`
      !! there. Taking the wrong number of tokens off the front shifts every value
      !! by one and reads the tensor into the centroid, so the count is per section
      !! rather than shared.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      type(efp_fragment_t), intent(inout) :: frag
      type(error_t), intent(inout) :: error

      integer :: start, finish, i, k, n_rec, at, stat
      character(len=:), allocatable :: joined, rest
      character(len=MAX_LINE) :: text
      real(dp) :: values(N_DYNAMIC_RECORD)

      start = 0
      do i = 1, n_lines
         text = adjustl(lines(i))
         ! Matched exactly: `DYNAMIC POLARIZABLE POINTS` and
         ! `LMOQQPOL DYNAMIC POLARIZABLE POINTS` both contain this name.
         if (trim(text) == "POLARIZABLE POINTS") then
            start = i + 1
            exit
         end if
      end do
      if (start == 0) return

      finish = start - 1
      do i = start, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         finish = i
      end do

      n_rec = 0
      i = start
      do while (i <= finish)
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) > 0) n_rec = n_rec + 1
      end do
      if (n_rec == 0) return

      frag%n_pol = n_rec
      allocate (frag%static_pol(3, 3, n_rec), frag%pol_points(3, n_rec))
      frag%static_pol = 0.0_dp
      i = start
      at = 0
      do while (i <= finish)
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         at = at + 1
         rest = strip_tokens(joined, 1)
         do k = 1, N_DYNAMIC_RECORD
            call next_number(rest, values(k), stat)
            if (stat /= 0) then
               call error%set(ERROR_VALIDATION, "efp: a polarizable point is short "// &
                              "of its twelve numbers")
               return
            end if
         end do
         frag%pol_points(:, at) = values(1:3)
         do k = 1, N_POL_SLOTS
            frag%static_pol(POL_ROW(k), POL_COL(k), at) = values(3 + k)
         end do
      end do
      frag%has_static_pol = .true.
   end subroutine read_static_pol

   subroutine read_dynamic(lines, n_lines, frag, error)
      !! `DYNAMIC POLARIZABLE POINTS`: a 3x3 tensor per localized orbital per
      !! frequency, which is what dispersion is built from
      !!
      !! The section's shape is a block per frequency, each block one record per
      !! localized orbital: a label, the orbital's centroid, then nine tensor
      !! components over continuation lines. Only the *first* record of a block
      !! carries the frequency, tagged on the end of its label line as
      !! `-- FOR W= 0.002792I A.U.`, so a new block is recognised by that tag
      !! appearing rather than by counting records.
      !!
      !! Absent is not an error: a potential written without the dynamic response
      !! is still usable for electrostatics.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      type(efp_fragment_t), intent(inout) :: frag
      type(error_t), intent(inout) :: error

      integer :: start, finish, i, k, n_rec, n_blocks, lmo, freq
      integer :: per_block
      character(len=:), allocatable :: joined, rest
      character(len=MAX_LINE) :: text
      real(dp) :: values(N_DYNAMIC_RECORD)
      integer :: stat
      logical :: new_block

      start = 0
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (index(text, "DYNAMIC POLARIZABLE POINTS") == 1) then
            start = i + 1
            exit
         end if
      end do
      if (start == 0) return

      finish = start - 1
      do i = start, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         finish = i
      end do

      ! First pass: how many records, and how many of them before the frequency tag
      ! reappears -- that count is the number of localized orbitals.
      n_rec = 0
      n_blocks = 0
      per_block = 0
      i = start
      do while (i <= finish)
         new_block = index(lines(i), "FOR W=") > 0
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         if (new_block) n_blocks = n_blocks + 1
         n_rec = n_rec + 1
         if (n_blocks == 1) per_block = per_block + 1
      end do
      if (n_rec == 0 .or. n_blocks == 0) return
      if (per_block*n_blocks /= n_rec) then
         call error%set(ERROR_VALIDATION, "efp: the dynamic polarizability section "// &
                        "does not hold the same number of orbitals at every frequency")
         return
      end if

      frag%n_lmo = per_block
      frag%n_freq = n_blocks
      allocate (frag%dyn_pol(3, 3, per_block, n_blocks))
      allocate (frag%centroids(3, per_block), frag%frequencies(n_blocks))
      frag%frequencies = 0.0_dp

      i = start
      freq = 0
      lmo = 0
      do while (i <= finish)
         new_block = index(lines(i), "FOR W=") > 0
         if (new_block) then
            freq = freq + 1
            lmo = 0
            call frequency_of(lines(i), frag%frequencies(freq))
         end if
         call join_dynamic(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         lmo = lmo + 1
         ! The label is two tokens, `CT` and its index, so three coordinates and
         ! nine tensor components follow from the third.
         rest = strip_tokens(joined, 2)
         do k = 1, N_DYNAMIC_RECORD
            call next_number(rest, values(k), stat)
            if (stat /= 0) then
               call error%set(ERROR_VALIDATION, "efp: a dynamic polarizability "// &
                              "record is short of its twelve numbers")
               return
            end if
         end do
         if (freq == 1) frag%centroids(:, lmo) = values(1:3)
         do k = 1, N_POL_SLOTS
            frag%dyn_pol(POL_ROW(k), POL_COL(k), lmo, freq) = values(3 + k)
         end do
      end do
      frag%has_dynamic = .true.
   end subroutine read_dynamic

   subroutine join_dynamic(lines, i, finish, joined)
      !! One dynamic-polarizability record: its label line plus its tensor lines
      !!
      !! **This section cannot use `join_record`.** There, a record is continued
      !! only where a line ends in `>`, and here the label line does not -- the
      !! tensor starts on the line *after* it, and the `>` marks continue within the
      !! tensor. Reading it the general way makes each tensor line a record of its
      !! own, which is how this announced itself: the record count stopped dividing
      !! by the frequency count.
      !!
      !! The frequency tag is dropped: everything from `--` onwards is a comment on
      !! the label line, not data.
      character(len=*), intent(in) :: lines(:)
      integer, intent(inout) :: i
      integer, intent(in) :: finish
      character(len=:), allocatable, intent(out) :: joined

      character(len=MAX_LINE) :: text
      integer :: cut
      logical :: continued

      joined = ""
      if (i > finish) return
      text = lines(i)
      i = i + 1
      cut = index(text, "--")
      if (cut > 0) then
         joined = trim(text(1:cut - 1))
      else
         joined = trim(text)
      end if

      do
         if (i > finish) exit
         text = lines(i)
         ! A new label line means this record is done; do not consume it.
         if (index(adjustl(text), "CT") == 1) exit
         i = i + 1
         cut = index(text, ">")
         continued = cut > 0
         if (continued) then
            joined = joined//" "//trim(text(1:cut - 1))
         else
            joined = joined//" "//trim(text)
         end if
         if (.not. continued) exit
      end do
      joined = adjustl(joined)
   end subroutine join_dynamic

   function strip_tokens(text, n) result(rest)
      !! `text` with its first `n` whitespace-separated tokens removed
      character(len=*), intent(in) :: text
      integer, intent(in) :: n
      character(len=:), allocatable :: rest

      integer :: k, space

      rest = adjustl(text)
      do k = 1, n
         space = index(trim(rest), " ")
         if (space == 0) then
            rest = ""
            return
         end if
         rest = adjustl(rest(space:))
      end do
   end function strip_tokens

   subroutine frequency_of(line, w)
      !! The imaginary frequency out of `-- FOR W= 0.002792I A.U.`
      character(len=*), intent(in) :: line
      real(dp), intent(out) :: w

      integer :: at, stop_at, stat

      w = 0.0_dp
      at = index(line, "FOR W=")
      if (at == 0) return
      stop_at = index(line(at:), "I")
      if (stop_at == 0) return
      read (line(at + 6:at + stop_at - 2), *, iostat=stat) w
      if (stat /= 0) w = 0.0_dp
   end subroutine frequency_of

   subroutine one_section(lines, n_lines, name, n_values, n_expected, values, error)
      !! A section that must be present and must carry one record per point
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      character(len=*), intent(in) :: name
      integer, intent(in) :: n_values, n_expected
      real(dp), allocatable, intent(out) :: values(:, :)
      type(error_t), intent(inout) :: error

      character(len=MAX_LINE), allocatable :: labels(:)
      integer :: n_rec

      call section_records(lines, n_lines, name, n_values, labels, values, n_rec, error)
      if (error%has_error()) return
      if (n_rec == 0) then
         call error%set(ERROR_VALIDATION, "efp: no "//name//" section")
         return
      end if
      call expect_points(n_rec, n_expected, name, error)
      if (allocated(labels)) deallocate (labels)
   end subroutine one_section

   subroutine expect_points(got, want, name, error)
      integer, intent(in) :: got, want
      character(len=*), intent(in) :: name
      type(error_t), intent(inout) :: error

      character(len=32) :: a, b

      if (got /= want) then
         write (a, "(I0)") got
         write (b, "(I0)") want
         call error%set(ERROR_VALIDATION, "efp: "//name//" carries "//trim(a)// &
                        " records but COORDINATES has "//trim(b)//" points")
      end if
   end subroutine expect_points

   subroutine section_records(lines, n_lines, name, n_values, labels, values, n_rec, error)
      !! Every record of one named section, continuations joined
      !!
      !! Returns `n_rec = 0` rather than an error when the section is absent, so an
      !! optional section costs the caller one branch.
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      character(len=*), intent(in) :: name
      integer, intent(in) :: n_values
      character(len=MAX_LINE), allocatable, intent(out) :: labels(:)
      real(dp), allocatable, intent(out) :: values(:, :)
      integer, intent(out) :: n_rec
      type(error_t), intent(inout) :: error

      integer :: i, start, finish, k, count_rec
      character(len=:), allocatable :: joined
      character(len=MAX_LINE) :: text

      n_rec = 0
      start = 0
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (trim(text) == trim(name) .or. index(trim(text), trim(name)//" ") == 1) then
            ! A header, not a record: records carry a label in column one after
            ! adjustment and this matched the section name outright.
            start = i + 1
            exit
         end if
      end do
      if (start == 0) then
         allocate (labels(0), values(n_values, 0))
         return
      end if

      finish = start - 1
      do i = start, n_lines
         if (trim(adjustl(lines(i))) == "STOP") exit
         finish = i
      end do

      ! First pass: count records, a record being a line that does not continue the
      ! one before it.
      count_rec = 0
      i = start
      do while (i <= finish)
         call join_record(lines, i, finish, joined)
         if (len_trim(joined) > 0) count_rec = count_rec + 1
      end do
      allocate (labels(count_rec), values(n_values, count_rec))
      values = 0.0_dp

      i = start
      k = 0
      do while (i <= finish)
         call join_record(lines, i, finish, joined)
         if (len_trim(joined) == 0) cycle
         k = k + 1
         call split_record(joined, n_values, labels(k), values(:, k), name, error)
         if (error%has_error()) return
      end do
      n_rec = k
   end subroutine section_records

   subroutine join_record(lines, i, finish, joined)
      !! One logical record starting at `i`, following `>` continuations
      !!
      !! `i` is advanced past everything consumed, so the caller's loop needs no
      !! separate bookkeeping for how many physical lines a record took.
      character(len=*), intent(in) :: lines(:)
      integer, intent(inout) :: i
      integer, intent(in) :: finish
      character(len=:), allocatable, intent(out) :: joined

      character(len=MAX_LINE) :: text
      integer :: cut

      joined = ""
      do
         if (i > finish) exit
         text = lines(i)
         i = i + 1
         cut = index(text, ">")
         if (cut > 0) then
            joined = joined//" "//trim(text(1:cut - 1))
         else
            joined = joined//" "//trim(text)
            exit
         end if
      end do
      joined = adjustl(joined)
   end subroutine join_record

   subroutine split_record(record, n_values, label, values, name, error)
      !! A record's label and its numbers
      character(len=*), intent(in) :: record
      integer, intent(in) :: n_values
      character(len=*), intent(out) :: label
      real(dp), intent(out) :: values(:)
      character(len=*), intent(in) :: name
      type(error_t), intent(inout) :: error

      integer :: space, stat, k
      character(len=:), allocatable :: rest
      character(len=32) :: a

      space = index(trim(record), " ")
      if (space == 0) then
         call error%set(ERROR_VALIDATION, "efp: "//name//": record '"// &
                        trim(record)//"' carries no values")
         return
      end if
      label = record(1:space - 1)
      rest = adjustl(record(space:))
      do k = 1, n_values
         call next_number(rest, values(k), stat)
         if (stat /= 0) then
            write (a, "(I0)") n_values
            call error%set(ERROR_VALIDATION, "efp: "//name//": record '"// &
                           trim(label)//"' needs "//trim(a)//" values")
            return
         end if
      end do
   end subroutine split_record

   subroutine next_number(text, value, stat)
      !! Read one number off the front of `text` and consume it
      character(len=:), allocatable, intent(inout) :: text
      real(dp), intent(out) :: value
      integer, intent(out) :: stat

      integer :: space

      value = 0.0_dp
      text = adjustl(text)
      if (len_trim(text) == 0) then
         stat = 1
         return
      end if
      space = index(trim(text), " ")
      if (space == 0) then
         read (text, *, iostat=stat) value
         text = ""
      else
         read (text(1:space - 1), *, iostat=stat) value
         text = adjustl(text(space:))
      end if
   end subroutine next_number

   subroutine fragment_name(lines, n_lines, name)
      !! `$FRAGNAME` without its `$`, or a placeholder if the file omits it
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      character(len=:), allocatable, intent(out) :: name

      integer :: i
      character(len=MAX_LINE) :: text

      name = "FRAGMENT"
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (text(1:1) == "$") then
            name = trim(text(2:))
            return
         end if
      end do
   end subroutine fragment_name

   subroutine multiplicity(lines, n_lines, mult)
      character(len=*), intent(in) :: lines(:)
      integer, intent(in) :: n_lines
      integer, intent(out) :: mult

      integer :: i, stat
      character(len=MAX_LINE) :: text

      mult = 1
      do i = 1, n_lines
         text = adjustl(lines(i))
         if (index(text, "MULTIPLICITY") == 1) then
            read (text(13:), *, iostat=stat) mult
            if (stat /= 0) mult = 1
            return
         end if
      end do
   end subroutine multiplicity

   subroutine slurp(path, lines, n_lines, error)
      !! The whole file as lines
      !!
      !! Read twice rather than grown: a potential is a few thousand lines, and one
      !! extra pass over it costs nothing against knowing the count up front.
      character(len=*), intent(in) :: path
      character(len=MAX_LINE), allocatable, intent(out) :: lines(:)
      integer, intent(out) :: n_lines
      type(error_t), intent(inout) :: error

      integer :: unit, stat, i
      character(len=MAX_LINE) :: text

      n_lines = 0
      open (newunit=unit, file=path, status="old", action="read", iostat=stat)
      if (stat /= 0) then
         call error%set(ERROR_IO, "efp: cannot open '"//trim(path)//"'")
         return
      end if
      do
         read (unit, "(A)", iostat=stat) text
         if (stat /= 0) exit
         n_lines = n_lines + 1
      end do
      rewind (unit)
      allocate (lines(max(n_lines, 1)))
      lines = ""
      do i = 1, n_lines
         read (unit, "(A)", iostat=stat) lines(i)
         if (stat /= 0) exit
      end do
      close (unit)
   end subroutine slurp

end module mqc_efp_read
