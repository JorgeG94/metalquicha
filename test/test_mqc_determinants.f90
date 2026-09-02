!! Occupation strings, addressing, and the single-excitation table
module test_mqc_determinants
   !! The first stage of a complete-active-space CI, and the one that can be
   !! checked completely without any physics at all. Everything asserted here
   !! follows from combinatorics and anticommutation: the counts are binomial
   !! coefficients, the addressing is a bijection, and the phases are the parity
   !! of how many electrons an operator has to step over.
   !!
   !! Two kinds of check, and both are needed. The *identities* -- round-trip
   !! addressing, row counts, phase antisymmetry -- are asserted over every
   !! string of several active spaces, so they are exhaustive rather than
   !! sampled, and they would catch an implementation that is self-consistently
   !! wrong in only one direction. The *reference values* are transcribed from
   !! `pyscf.fci.cistring` and catch the case the identities cannot: a
   !! convention that is internally consistent and different from everyone
   !! else's, which would go unnoticed here and produce a silently wrong CI
   !! energy three stages later.
   !!
   !! PySCF indexes orbitals and addresses from zero and this code from one, so
   !! every transcribed index below is its value plus one. The phases are not
   !! shifted -- a sign is a sign.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: int64
   use mqc_error, only: error_t
   use mqc_determinants, only: n_strings, generate_strings, string_address, &
                               address_to_string, excitation_phase, &
                               link_table_t, build_link_table, MAX_ACTIVE_ORBITALS
   implicit none
   private

   public :: collect_mqc_determinants_tests

contains

   subroutine collect_mqc_determinants_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("string_counts", test_counts), &
                  new_unittest("strings_match_pyscf", test_strings), &
                  new_unittest("addressing_round_trips", test_round_trip), &
                  new_unittest("excitation_phases", test_phases), &
                  new_unittest("link_table_against_pyscf", test_link_reference), &
                  new_unittest("link_table_identities", test_link_identities), &
                  new_unittest("refusals", test_refusals) &
                  ]
   end subroutine collect_mqc_determinants_tests

   subroutine test_counts(error)
      !! The number of strings is a binomial coefficient
      type(error_type), allocatable, intent(out) :: error

      call check(error, n_strings(4, 2) == 6_int64, "C(4,2) = 6")
      if (allocated(error)) return
      call check(error, n_strings(6, 3) == 20_int64, "C(6,3) = 20")
      if (allocated(error)) return
      call check(error, n_strings(1, 0) == 1_int64, "there is one empty string")
      if (allocated(error)) return
      call check(error, n_strings(8, 8) == 1_int64, "a full shell has one string")
      if (allocated(error)) return
      call check(error, n_strings(4, 5) == 0_int64, &
                 "more electrons than orbitals is no strings, not an error here")
      if (allocated(error)) return

      ! The determinant counts the active spaces in the literature are quoted
      ! by. CAS(12,12) is 924 strings per spin and 853776 determinants, which is
      ! the size the whole design has to survive.
      call check(error, n_strings(10, 5) == 252_int64, "CAS(10,10) has 252 strings")
      if (allocated(error)) return
      call check(error, n_strings(12, 6) == 924_int64, "CAS(12,12) has 924 strings")
      if (allocated(error)) return
      call check(error, n_strings(14, 7) == 3432_int64, "CAS(14,14) has 3432 strings")
      if (allocated(error)) return

      ! Computed by multiply-and-divide rather than from factorials, which is
      ! the only reason this does not overflow: 40! is about 8e47 and the answer
      ! fits in a 64-bit integer with room to spare.
      call check(error, n_strings(40, 20) == 137846528820_int64, &
                 "C(40,20) should come out exactly, where a factorial would have "// &
                 "overflowed long before")
   end subroutine test_counts

   subroutine test_strings(error)
      !! The strings themselves, against PySCF
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer(int64), allocatable :: strings(:)
      integer :: i

      !! `cistring.make_strings(range(4), 2)`
      integer(int64), parameter :: FOUR_TWO(6) = &
                                   [3_int64, 5_int64, 6_int64, 9_int64, 10_int64, 12_int64]
      !! `cistring.make_strings(range(6), 3)`
      integer(int64), parameter :: SIX_THREE(20) = &
                                   [7_int64, 11_int64, 13_int64, 14_int64, 19_int64, &
                                    21_int64, 22_int64, 25_int64, 26_int64, 28_int64, &
                                    35_int64, 37_int64, 38_int64, 41_int64, 42_int64, &
                                    44_int64, 49_int64, 50_int64, 52_int64, 56_int64]

      call generate_strings(4, 2, strings, err)
      call check(error,.not. err%has_error(), "two electrons in four orbitals")
      if (allocated(error)) return
      call check(error, size(strings), 6, "six strings")
      if (allocated(error)) return
      do i = 1, 6
         call check(error, strings(i) == FOUR_TWO(i), "string "//char(48 + i))
         if (allocated(error)) return
      end do
      deallocate (strings)

      call generate_strings(6, 3, strings, err)
      call check(error,.not. err%has_error(), "three electrons in six orbitals")
      if (allocated(error)) return
      call check(error, size(strings), 20, "twenty strings")
      if (allocated(error)) return
      call check(error, all(strings == SIX_THREE), &
                 "every string of CAS(6,6) alpha space, against PySCF")
      if (allocated(error)) return

      ! Strictly ascending, which is what makes the generator's order and the
      ! addressing agree. Asserted separately because Gosper's hack producing
      ! something out of order is a different failure from it producing the
      ! wrong set.
      do i = 2, size(strings)
         call check(error, strings(i) > strings(i - 1), &
                    "strings should come out strictly ascending")
         if (allocated(error)) return
      end do
   end subroutine test_strings

   subroutine test_round_trip(error)
      !! Address and string are inverse, for every string of several spaces
      !!
      !! Exhaustive rather than sampled: over 3400 strings across the six spaces
      !! below, each checked in both directions. The address is what indexes the
      !! CI vector, so a collision or a gap here is not a small error, it is a
      !! different Hilbert space.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      integer(int64), allocatable :: strings(:)
      integer :: space, i, norb, nelec
      integer, parameter :: ORBITALS(6) = [4, 6, 8, 8, 10, 14]
      integer, parameter :: ELECTRONS(6) = [2, 3, 1, 4, 5, 7]

      do space = 1, 6
         norb = ORBITALS(space)
         nelec = ELECTRONS(space)
         call generate_strings(norb, nelec, strings, err)
         call check(error,.not. err%has_error(), "the space should build")
         if (allocated(error)) return

         do i = 1, size(strings)
            call check(error, string_address(norb, nelec, strings(i)), i, &
                       "the address of the i'th string should be i")
            if (allocated(error)) return
            call check(error, address_to_string(norb, nelec, i) == strings(i), &
                       "and the string at address i should be the i'th string")
            if (allocated(error)) return
            call check(error, int(popcnt(strings(i))), nelec, &
                       "every string should hold exactly the electron count")
            if (allocated(error)) return
         end do
         deallocate (strings)
      end do
   end subroutine test_round_trip

   subroutine test_phases(error)
      !! The sign left behind by a single excitation
      type(error_type), allocatable, intent(out) :: error
      integer(int64) :: string

      ! Transcribed from `cistring.cre_des_sign`, orbitals shifted by one.
      ! p=2 q=0 on 0b0011: one electron sits between them, so odd parity.
      call check(error, excitation_phase(3, 1, 3_int64), -1, &
                 "stepping over one electron changes the sign")
      if (allocated(error)) return
      ! p=3 q=0 on 0b0111: two electrons between, even parity.
      call check(error, excitation_phase(4, 1, 7_int64), 1, &
                 "stepping over two electrons does not")
      if (allocated(error)) return
      ! Blocked both ways.
      call check(error, excitation_phase(4, 2, 11_int64), 0, &
                 "creating into an occupied orbital annihilates the string")
      if (allocated(error)) return
      call check(error, excitation_phase(1, 3, 5_int64), 0, &
                 "so does annihilating an empty one")
      if (allocated(error)) return

      ! The number operator, where this and PySCF deliberately differ: PySCF
      ! returns 1 without looking, because it only ever asks about occupied
      ! orbitals. Asked about an empty one, the answer is zero.
      call check(error, excitation_phase(2, 2, 3_int64), 1, &
                 "the number operator on an occupied orbital")
      if (allocated(error)) return
      call check(error, excitation_phase(3, 3, 3_int64), 0, &
                 "and on an empty one")
      if (allocated(error)) return

      ! Adjacent orbitals have nothing between them whichever way round they
      ! are, and direction cannot matter to the count.
      string = 5_int64      ! orbitals 1 and 3 occupied
      call check(error, excitation_phase(2, 1, string), 1, &
                 "adjacent orbitals, upward")
      if (allocated(error)) return
      call check(error, excitation_phase(2, 3, string), 1, &
                 "adjacent orbitals, downward")
   end subroutine test_phases

   subroutine test_link_reference(error)
      !! The whole excitation table for CAS(4,2), against PySCF
      !!
      !! Small enough to write out completely, which is the point: every row of
      !! every string, so there is nowhere for a convention difference to hide.
      !! `gen_linkstr_index(range(4), 2)` gives a (6, 6, 4) array and this is it,
      !! with orbitals and addresses shifted by one.
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: table
      integer :: istr, row

      !! (cre, des, dest, phase) for each of 6 rows of each of 6 strings.
      integer, parameter :: REF(4, 6, 6) = reshape([ &
                                                   ! str 1, 0b0011
                                               1, 1, 1, 1, 2, 2, 1, 1, 3, 1, 3, -1, 4, 1, 5, -1, 3, 2, 2, 1, 4, 2, 4, 1, &
                                                   ! str 2, 0b0101
                                                1, 1, 2, 1, 3, 3, 2, 1, 2, 1, 3, 1, 4, 1, 6, -1, 2, 3, 1, 1, 4, 3, 4, 1, &
                                                   ! str 3, 0b0110
                                               2, 2, 3, 1, 3, 3, 3, 1, 1, 2, 2, 1, 4, 2, 6, -1, 1, 3, 1, -1, 4, 3, 5, 1, &
                                                   ! str 4, 0b1001
                                                 1, 1, 4, 1, 4, 4, 4, 1, 2, 1, 5, 1, 3, 1, 6, 1, 2, 4, 1, 1, 3, 4, 2, 1, &
                                                   ! str 5, 0b1010
                                                2, 2, 5, 1, 4, 4, 5, 1, 1, 2, 4, 1, 3, 2, 6, 1, 1, 4, 1, -1, 3, 4, 3, 1, &
                                                   ! str 6, 0b1100
                                              3, 3, 6, 1, 4, 4, 6, 1, 1, 3, 4, 1, 2, 3, 5, 1, 1, 4, 2, -1, 2, 4, 3, -1], &
                                                   [4, 6, 6])

      call build_link_table(4, 2, table, err)
      call check(error,.not. err%has_error(), "the table should build")
      if (allocated(error)) return
      call check(error, table%n_strings, 6, "six strings")
      if (allocated(error)) return
      call check(error, table%n_rows, 6, "two diagonal rows and four excitations")
      if (allocated(error)) return

      do istr = 1, 6
         do row = 1, 6
            call check(error, table%cre(row, istr), REF(1, row, istr), "cre")
            if (allocated(error)) return
            call check(error, table%des(row, istr), REF(2, row, istr), "des")
            if (allocated(error)) return
            call check(error, table%dest(row, istr), REF(3, row, istr), "dest")
            if (allocated(error)) return
            call check(error, table%phase(row, istr), REF(4, row, istr), "phase")
            if (allocated(error)) return
         end do
      end do

      call table%destroy()
   end subroutine test_link_reference

   subroutine test_link_identities(error)
      !! Properties the table has to have, over several active spaces
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: table
      integer(int64), allocatable :: strings(:)
      integer(int64) :: source, arrived
      integer :: space, istr, row, norb, nelec, diagonal
      integer, parameter :: ORBITALS(4) = [4, 6, 8, 10]
      integer, parameter :: ELECTRONS(4) = [2, 3, 4, 5]

      do space = 1, 4
         norb = ORBITALS(space)
         nelec = ELECTRONS(space)
         call build_link_table(norb, nelec, table, err)
         call check(error,.not. err%has_error(), "the table should build")
         if (allocated(error)) return
         call generate_strings(norb, nelec, strings, err)
         if (allocated(error)) return

         call check(error, table%n_rows, nelec*(1 + norb - nelec), &
                    "one row per occupied orbital plus one per occupied-virtual pair")
         if (allocated(error)) return

         do istr = 1, table%n_strings
            source = strings(istr)
            diagonal = 0
            do row = 1, table%n_rows
               ! No row may be a phase of zero: every entry describes an
               ! excitation that actually connects two strings, and a zero here
               ! would mean the builder tabulated one that does not exist.
               call check(error, abs(table%phase(row, istr)), 1, &
                          "every tabulated excitation should be a real one")
               if (allocated(error)) return

               ! The destination is what the operators actually produce.
               arrived = ibset(ibclr(source, table%des(row, istr) - 1), &
                               table%cre(row, istr) - 1)
               call check(error, string_address(norb, nelec, arrived), &
                          table%dest(row, istr), &
                          "the tabulated address should be the address of the "// &
                          "string the excitation produces")
               if (allocated(error)) return

               ! The orbital annihilated from must be occupied in the source,
               ! and the one created into must not be, unless they are the same.
               call check(error, btest(source, table%des(row, istr) - 1), &
                          "the annihilated orbital should be occupied")
               if (allocated(error)) return
               if (table%cre(row, istr) == table%des(row, istr)) then
                  diagonal = diagonal + 1
                  call check(error, table%dest(row, istr), istr, &
                             "a number operator maps a string to itself")
                  if (allocated(error)) return
               else
                  call check(error,.not. btest(source, table%cre(row, istr) - 1), &
                             "the created orbital should be empty")
                  if (allocated(error)) return
               end if
            end do
            call check(error, diagonal, nelec, &
                       "there should be exactly one number operator per electron")
            if (allocated(error)) return
         end do

         call table%destroy()
         deallocate (strings)
      end do
   end subroutine test_link_identities

   subroutine test_refusals(error)
      !! Spaces this cannot address
      type(error_type), allocatable, intent(out) :: error
      type(error_t) :: err
      type(link_table_t) :: table
      integer(int64), allocatable :: strings(:)

      ! More electrons of one spin than there are orbitals to hold them. Worth
      ! refusing loudly because the usual cause is passing the total electron
      ! count where the per-spin one belongs, which for a closed shell is a
      ! factor of two and produces a smaller space rather than an obviously
      ! broken one.
      call generate_strings(4, 6, strings, err)
      call check(error, err%has_error(), &
                 "six electrons of one spin in four orbitals should be refused")
      if (allocated(error)) return
      call err%clear()

      ! Past the width of the representation.
      call build_link_table(MAX_ACTIVE_ORBITALS + 1, 2, table, err)
      call check(error, err%has_error(), &
                 "an active space wider than a 64-bit string should be refused")
      if (allocated(error)) return
      call err%clear()

      ! Addressable in principle, unstorable in practice: C(60,30) overflows a
      ! default integer address long before memory becomes the objection.
      call build_link_table(60, 30, table, err)
      call check(error, err%has_error(), &
                 "a space with more strings than a default integer can address "// &
                 "should be refused rather than silently wrapping")
      call err%clear()
   end subroutine test_refusals

end module test_mqc_determinants

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_determinants, only: collect_mqc_determinants_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_determinants", collect_mqc_determinants_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester
