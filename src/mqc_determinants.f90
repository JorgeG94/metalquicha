!! Occupation strings, their addresses, and the excitations connecting them
module mqc_determinants
   !! The bookkeeping a complete-active-space CI is built on, and nothing else:
   !! no integrals, no orbitals, no energies. A determinant in a CAS is a pair of
   !! *strings*, one per spin, each recording which of the active orbitals that
   !! spin occupies. This module generates them, numbers them, and tabulates the
   !! single excitations between them.
   !!
   !! **Not the strings in `pic_strings`**, which manipulate text. The word is
   !! overloaded and this is the other sense: a bit pattern over active orbitals,
   !! bit `k` set meaning orbital `k + 1` is occupied.
   !!
   !! Everything here is exactly checkable without any physics -- the counts are
   !! binomial coefficients, the addressing is a bijection that has to round
   !! trip, and the excitation phases follow from anticommutation alone. That is
   !! the reason to build it first and on its own. A wrong phase in the table
   !! below does not announce itself later; it produces a CI energy that is
   !! merely wrong, with nothing obviously out of place, and by then there are
   !! three other stages to suspect.
   !!
   !! The conventions -- string order, address order, and the layout of the
   !! excitation table -- follow PySCF's `pyscf/fci/cistring.py`, deliberately.
   !! They are as good as any other choice and they are the ones the reference
   !! this will be validated against uses, which means intermediate quantities
   !! can be compared and not just final energies.
   !!
   !! Orbital indices and addresses are 1-based here, as everywhere else in this
   !! code, where PySCF's are 0-based. Bit `k` of a string is orbital `k + 1`.
   use pic_types, only: int64, int32
   use pic_io, only: to_char
   use mqc_error, only: error_t, ERROR_VALIDATION
   use mqc_math_utils, only: binomial
   implicit none
   private

   public :: MAX_ACTIVE_ORBITALS
   public :: n_strings
   public :: generate_strings
   public :: string_address
   public :: address_to_string
   public :: excitation_phase
   public :: link_table_t
   public :: build_link_table

   integer, parameter :: MAX_ACTIVE_ORBITALS = 63
      !! A string is one `int64`, so this is where the representation stops.
      !! Nothing comes close: the determinant count at CAS(63,31) exceeds the
      !! number of atoms in the observable universe, and `build_link_table`
      !! refuses far below it.

   type :: link_table_t
      !! Every single excitation out of every string
      !!
      !! Column `s` lists the excitations `E_ai` acting on string `s`:
      !! `cre(r, s)` and `des(r, s)` are the orbitals created in and annihilated
      !! from, `dest(r, s)` is the address of the string that results, and
      !! `phase(r, s)` is the sign anticommutation leaves behind. So
      !!
      !!     a^dagger_{cre} a_{des} |s> = phase * |dest>
      !!
      !! The first `n_electrons` rows are the diagonal ones, `cre == des` over
      !! the occupied orbitals, which map a string to itself with phase `+1`.
      !! They are in the table rather than special-cased because the number
      !! operator terms of the Hamiltonian need them and a caller that had to
      !! remember to add them separately would eventually not.
      !!
      !! Stored with the row index first so one string's excitations are
      !! contiguous, which is the order a sigma build walks them in.
      integer :: n_orbitals = 0
      integer :: n_electrons = 0
      integer :: n_strings = 0
      integer :: n_rows = 0
         !! `n_electrons * (1 + n_orbitals - n_electrons)` -- the diagonal
         !! entries plus one per occupied-virtual pair.
      integer, allocatable :: cre(:, :)     !! (n_rows, n_strings), 1-based orbital
      integer, allocatable :: des(:, :)     !! (n_rows, n_strings), 1-based orbital
      integer, allocatable :: dest(:, :)    !! (n_rows, n_strings), 1-based address
      integer, allocatable :: phase(:, :)   !! (n_rows, n_strings), +1 or -1
   contains
      procedure :: destroy => link_table_destroy
   end type link_table_t

contains

   pure function n_strings(n_orbitals, n_electrons) result(count)
      !! How many strings put `n_electrons` into `n_orbitals`
      integer, intent(in) :: n_orbitals, n_electrons
      integer(int64) :: count

      count = binomial(n_orbitals, n_electrons)
   end function n_strings

   subroutine check_space(n_orbitals, n_electrons, count, error)
      !! The three ways an active space can be one this module cannot address
      integer, intent(in) :: n_orbitals, n_electrons
      integer(int64), intent(out) :: count
      type(error_t), intent(inout) :: error

      count = 0_int64
      if (error%has_error()) return

      if (n_electrons < 0 .or. n_orbitals < 0 .or. n_electrons > n_orbitals) then
         call error%set(ERROR_VALIDATION, "an active space of "// &
                        to_char(n_electrons)//" electrons of one spin in "// &
                        to_char(n_orbitals)//" orbitals is not a space. The count "// &
                        "here is per spin, so it is half the electrons, not all of them.")
         return
      end if
      if (n_orbitals > MAX_ACTIVE_ORBITALS) then
         call error%set(ERROR_VALIDATION, "the active space has "// &
                        to_char(n_orbitals)//" orbitals and a string is one 64-bit "// &
                        "integer, so "//to_char(MAX_ACTIVE_ORBITALS)//" is the limit "// &
                        "of the representation. Nothing of that size is computable "// &
                        "by this route anyway.")
         return
      end if

      count = binomial(n_orbitals, n_electrons)
      if (count > int(huge(1_int32), int64)) then
         call error%set(ERROR_VALIDATION, "this active space has more strings than "// &
                        "a default integer can address. It is far past what the "// &
                        "determinant expansion could be stored for in any case.")
         count = 0_int64
      end if
   end subroutine check_space

   subroutine generate_strings(n_orbitals, n_electrons, strings, error)
      !! Every string, in address order
      !!
      !! Address order is ascending numerical order of the bit pattern, so the
      !! generator is Gosper's: given one integer, produce the next larger with
      !! the same number of bits set. That the two orders coincide is not a
      !! coincidence to rely on silently -- `address_round_trip` in the tests
      !! asserts it for every string of several spaces.
      integer, intent(in) :: n_orbitals, n_electrons
      integer(int64), allocatable, intent(out) :: strings(:)
      type(error_t), intent(inout) :: error

      integer(int64) :: count, s, lowest, ripple, ones
      integer :: i

      call check_space(n_orbitals, n_electrons, count, error)
      if (error%has_error()) return

      allocate (strings(count))
      if (count == 0) return

      ! The lowest string fills the lowest orbitals.
      s = shiftl(1_int64, n_electrons) - 1_int64
      do i = 1, int(count)
         strings(i) = s
         if (i == int(count)) exit
         ! Gosper's hack. `lowest` is the lowest set bit, `ripple` carries it
         ! upward, and the shifted difference puts the bits that were displaced
         ! back at the bottom.
         lowest = iand(s, -s)
         ripple = s + lowest
         ones = shiftr(ieor(s, ripple), 2)/lowest
         s = ior(ripple, ones)
      end do
   end subroutine generate_strings

   pure function string_address(n_orbitals, n_electrons, string) result(address)
      !! Where a string sits in the list, 1-based
      !!
      !! Walk the orbitals from the top down. Each occupied orbital `k` found
      !! with `m` electrons still to place contributes `C(k, m)` -- the number
      !! of strings that agree with this one above `k` and put a lower orbital
      !! there instead, all of which sort first.
      !!
      !! No error argument, deliberately: this is called once per row of the
      !! excitation table and the caller has already established the space is
      !! valid. A string with the wrong number of bits set returns a meaningless
      !! address rather than complaining.
      integer, intent(in) :: n_orbitals, n_electrons
      integer(int64), intent(in) :: string
      integer :: address

      integer(int64) :: offset
      integer :: orbital, remaining

      offset = 0_int64
      remaining = n_electrons
      do orbital = n_orbitals - 1, 0, -1
         if (remaining == 0) exit
         if (btest(string, orbital)) then
            offset = offset + binomial(orbital, remaining)
            remaining = remaining - 1
         end if
      end do
      address = int(offset) + 1
   end function string_address

   pure function address_to_string(n_orbitals, n_electrons, address) result(string)
      !! The inverse of `string_address`
      integer, intent(in) :: n_orbitals, n_electrons
      integer, intent(in) :: address
      integer(int64) :: string

      integer(int64) :: remainder, cumulative
      integer :: orbitals, electrons, orbital

      string = 0_int64
      remainder = int(address - 1, int64)
      orbitals = n_orbitals
      electrons = n_electrons

      do
         ! Nothing left to choose: the remaining electrons fill the lowest
         ! orbitals, which is the first string of the sub-problem.
         if (remainder == 0_int64 .or. electrons == 0 .or. electrons == orbitals) then
            if (electrons > 0) string = ior(string, shiftl(1_int64, electrons) - 1_int64)
            return
         end if
         do orbital = orbitals - 1, 0, -1
            cumulative = binomial(orbital, electrons)
            if (cumulative <= remainder) then
               string = ibset(string, orbital)
               remainder = remainder - cumulative
               electrons = electrons - 1
               orbitals = orbital
               exit
            end if
         end do
      end do
   end function address_to_string

   pure function excitation_phase(create, destroy, string) result(phase)
      !! The sign of `a^dagger_create a_destroy |string>`, or zero
      !!
      !! Zero when the excitation annihilates the string: creating into an
      !! orbital already occupied, or annihilating one that is empty. Otherwise
      !! the sign is `(-1)` to the number of electrons sitting strictly between
      !! the two orbitals, which is how many anticommutations it takes to move
      !! the operator into place.
      !!
      !! Orbitals are 1-based. `create == destroy` is the number operator, which
      !! is `+1` on an occupied orbital and zero on an empty one -- PySCF's
      !! `cre_des_sign` returns `1` unconditionally there because it is only
      !! ever asked about occupied orbitals, and this returns the answer for the
      !! question as asked instead.
      integer, intent(in) :: create, destroy
      integer(int64), intent(in) :: string
      integer :: phase

      integer(int64) :: between
      integer :: low, high

      phase = 0
      if (create == destroy) then
         if (btest(string, create - 1)) phase = 1
         return
      end if
      if (btest(string, create - 1)) return
      if (.not. btest(string, destroy - 1)) return

      low = min(create, destroy)
      high = max(create, destroy)
      ! Bits strictly between the two orbitals: set the ones below `high`, clear
      ! the ones at or below `low`.
      between = shiftl(1_int64, high - 1) - shiftl(1_int64, low)
      phase = 1
      if (mod(popcnt(iand(string, between)), 2) == 1) phase = -1
   end function excitation_phase

   subroutine build_link_table(n_orbitals, n_electrons, table, error)
      !! Tabulate every single excitation out of every string
      !!
      !! Row order follows PySCF: the diagonal entries first, then one block per
      !! occupied orbital in ascending order, each listing the virtual orbitals
      !! in ascending order. The order carries no meaning -- a contraction sums
      !! over rows -- but matching it makes the table directly comparable
      !! against the reference implementation, which is worth more than the
      !! freedom to choose.
      integer, intent(in) :: n_orbitals, n_electrons
      type(link_table_t), intent(out) :: table
      type(error_t), intent(inout) :: error

      integer(int64), allocatable :: strings(:)
      integer(int64) :: count, string, excited
      integer :: n_virtual, istr, row, occupied, virtual

      call check_space(n_orbitals, n_electrons, count, error)
      if (error%has_error()) return

      call generate_strings(n_orbitals, n_electrons, strings, error)
      if (error%has_error()) return

      n_virtual = n_orbitals - n_electrons
      table%n_orbitals = n_orbitals
      table%n_electrons = n_electrons
      table%n_strings = int(count)
      table%n_rows = n_electrons*(1 + n_virtual)

      allocate (table%cre(table%n_rows, table%n_strings))
      allocate (table%des(table%n_rows, table%n_strings))
      allocate (table%dest(table%n_rows, table%n_strings))
      allocate (table%phase(table%n_rows, table%n_strings))

      do istr = 1, table%n_strings
         string = strings(istr)

         ! The number operators: every occupied orbital, mapping the string to
         ! itself.
         row = 0
         do occupied = 1, n_orbitals
            if (.not. btest(string, occupied - 1)) cycle
            row = row + 1
            table%cre(row, istr) = occupied
            table%des(row, istr) = occupied
            table%dest(row, istr) = istr
            table%phase(row, istr) = 1
         end do

         ! Then every occupied-to-virtual single excitation.
         do occupied = 1, n_orbitals
            if (.not. btest(string, occupied - 1)) cycle
            do virtual = 1, n_orbitals
               if (btest(string, virtual - 1)) cycle
               row = row + 1
               excited = ibset(ibclr(string, occupied - 1), virtual - 1)
               table%cre(row, istr) = virtual
               table%des(row, istr) = occupied
               table%dest(row, istr) = string_address(n_orbitals, n_electrons, excited)
               table%phase(row, istr) = excitation_phase(virtual, occupied, string)
            end do
         end do
      end do

      deallocate (strings)
   end subroutine build_link_table

   subroutine link_table_destroy(this)
      class(link_table_t), intent(inout) :: this

      if (allocated(this%cre)) deallocate (this%cre)
      if (allocated(this%des)) deallocate (this%des)
      if (allocated(this%dest)) deallocate (this%dest)
      if (allocated(this%phase)) deallocate (this%phase)
      this%n_orbitals = 0
      this%n_electrons = 0
      this%n_strings = 0
      this%n_rows = 0
   end subroutine link_table_destroy

end module mqc_determinants
