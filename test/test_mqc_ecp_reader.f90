module test_mqc_ecp_reader
   !! Tests the BSE JSON effective-core-potential reader.
   !!
   !! The values below are read off `basis_sets/def2-ecp.json` by hand. They
   !! are worth stating explicitly rather than derived, because the two things
   !! most likely to go wrong here are silent: picking the wrong channel as the
   !! local one, and losing digits from the coefficients. Neither fails, both
   !! give a wrong potential.
   !!
   !! Rubidium is the lightest element def2-ECP covers, and gold is a case
   !! where the core is large enough that getting `ecp_electrons` wrong would
   !! be unmistakable in the electron count.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_ecp, only: atomic_ecp_type, molecular_ecp_type
   use mqc_json_ecp_reader, only: read_json_ecp_element, build_molecular_ecp_json
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_ecp_reader_tests

   character(len=*), parameter :: ECP_FILE = "../basis_sets/def2-ecp.json"
      !! Tests run with their source directory as the working directory

contains

   subroutine collect_mqc_ecp_reader_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("rubidium_channels", test_rubidium), &
                  new_unittest("local_channel_is_highest_l", test_local_is_highest), &
                  new_unittest("gold_core_electrons", test_gold), &
                  new_unittest("light_element_has_no_ecp", test_light_element), &
                  new_unittest("molecular_set_is_per_atom", test_molecular), &
                  new_unittest("core_electrons_sum_over_atoms", test_core_sum), &
                  new_unittest("missing_file_is_an_error", test_missing_file) &
                  ]
   end subroutine collect_mqc_ecp_reader_tests

   subroutine test_rubidium(error)
      !! Channel count, angular momenta and the first primitive of each
      type(error_type), allocatable, intent(out) :: error
      type(atomic_ecp_type) :: ecp
      type(error_t) :: read_error

      call read_json_ecp_element(ECP_FILE, "Rb", ecp, read_error)
      call check(error,.not. read_error%has_error(), read_error%get_message())
      if (allocated(error)) return

      call check(error, ecp%has_ecp, "Rb should carry an ECP")
      if (allocated(error)) return
      call check(error, ecp%core_electrons, 28, "Rb ECP replaces 28 electrons")
      if (allocated(error)) return

      ! def2-ECP gives Rb four channels: l = 3 local, plus 0, 1, 2 projected.
      call check(error, ecp%n_projected, 3, "Rb should have three projected channels")
      if (allocated(error)) return
      call check(error, ecp%local%ang_mom, 3, "Rb local channel is l = 3")
      if (allocated(error)) return

      ! The single l=3 primitive, read straight out of the file.
      call check(error, ecp%local%nprim, 1, "Rb local channel has one primitive")
      if (allocated(error)) return
      call check(error, ecp%local%radial_powers(1), 2, "Rb local r exponent is 2")
      if (allocated(error)) return
      call check(error, close_enough(ecp%local%exponents(1), 3.8431140_dp), &
                 "Rb local exponent should be 3.8431140")
      if (allocated(error)) return
      call check(error, close_enough(ecp%local%coefficients(1), -12.3169000_dp), &
                 "Rb local coefficient should be -12.3169000")
      if (allocated(error)) return

      ! The l=0 channel has three primitives; check the widest-magnitude one so
      ! a truncated read shows up.
      call check(error, projected_of(ecp, 0) > 0, "Rb should have an l = 0 channel")
      if (allocated(error)) return
      call check(error, ecp%projected(projected_of(ecp, 0))%nprim, 3, &
                 "Rb l = 0 channel has three primitives")
      if (allocated(error)) return
      call check(error, close_enough(ecp%projected(projected_of(ecp, 0))%coefficients(1), &
                                     89.5001980_dp), &
                 "Rb l = 0 leading coefficient should be 89.5001980")

      call ecp%destroy()
   end subroutine test_rubidium

   subroutine test_local_is_highest(error)
      !! The local channel must carry a higher l than any projected one
      !!
      !! This is the assumption the whole semi-local form rests on. Selecting
      !! the wrong channel produces a potential that is wrong everywhere and
      !! complains nowhere.
      type(error_type), allocatable, intent(out) :: error
      type(atomic_ecp_type) :: ecp
      type(error_t) :: read_error
      character(len=2) :: symbols(4)
      integer :: isym, iproj

      symbols = ["Rb", "Ag", "Au", "Rn"]
      do isym = 1, size(symbols)
         call read_json_ecp_element(ECP_FILE, trim(symbols(isym)), ecp, read_error)
         call check(error,.not. read_error%has_error(), read_error%get_message())
         if (allocated(error)) return
         call check(error, ecp%has_ecp, trim(symbols(isym))//" should carry an ECP")
         if (allocated(error)) return

         do iproj = 1, ecp%n_projected
            call check(error, ecp%projected(iproj)%ang_mom < ecp%local%ang_mom, &
                       trim(symbols(isym))//": a projected channel has l >= the local one")
            if (allocated(error)) return
         end do
         call ecp%destroy()
      end do
   end subroutine test_local_is_highest

   subroutine test_gold(error)
      !! A large core, where an error in the count would be obvious
      type(error_type), allocatable, intent(out) :: error
      type(atomic_ecp_type) :: ecp
      type(error_t) :: read_error

      call read_json_ecp_element(ECP_FILE, "Au", ecp, read_error)
      call check(error,.not. read_error%has_error(), read_error%get_message())
      if (allocated(error)) return
      call check(error, ecp%core_electrons, 60, "Au ECP replaces 60 electrons")
      if (allocated(error)) return
      ! Z = 79 less a 60-electron core leaves 19 valence electrons.
      call check(error, 79 - ecp%core_electrons, 19, &
                 "Au should be left with 19 valence electrons")

      call ecp%destroy()
   end subroutine test_gold

   subroutine test_light_element(error)
      !! An element the file does not cover is not an error
      type(error_type), allocatable, intent(out) :: error
      type(atomic_ecp_type) :: ecp
      type(error_t) :: read_error

      call read_json_ecp_element(ECP_FILE, "O", ecp, read_error)
      call check(error,.not. read_error%has_error(), &
                 "an element without an ECP should not be an error: "// &
                 read_error%get_message())
      if (allocated(error)) return
      call check(error,.not. ecp%has_ecp, "oxygen should carry no ECP")
      if (allocated(error)) return
      call check(error, ecp%core_electrons, 0, "an absent ECP replaces no electrons")

      call ecp%destroy()
   end subroutine test_light_element

   subroutine test_molecular(error)
      !! One entry per atom, in atom order, mixing ECP and non-ECP elements
      type(error_type), allocatable, intent(out) :: error
      type(molecular_ecp_type) :: mol_ecp
      type(error_t) :: read_error
      character(len=2) :: symbols(4)

      ! Au-O-Au-H: heavy, light, heavy again, light.
      symbols = ["Au", "O ", "Au", "H "]
      call build_molecular_ecp_json(ECP_FILE, symbols, mol_ecp, read_error)
      call check(error,.not. read_error%has_error(), read_error%get_message())
      if (allocated(error)) return

      call check(error, mol_ecp%natoms, 4, "one entry per atom")
      if (allocated(error)) return
      call check(error, mol_ecp%n_ecp_atoms, 2, "two of the four atoms carry an ECP")
      if (allocated(error)) return

      ! Entries must line up with atom order, not element order -- everything
      ! downstream indexes coordinates and charges the same way.
      call check(error, mol_ecp%atoms(1)%has_ecp, "atom 1 (Au) should carry an ECP")
      if (allocated(error)) return
      call check(error,.not. mol_ecp%atoms(2)%has_ecp, "atom 2 (O) should not")
      if (allocated(error)) return
      call check(error, mol_ecp%atoms(3)%has_ecp, "atom 3 (Au) should carry an ECP")
      if (allocated(error)) return
      call check(error,.not. mol_ecp%atoms(4)%has_ecp, "atom 4 (H) should not")
      if (allocated(error)) return

      ! The repeated element must be a real copy, not a shared or empty one.
      call check(error, mol_ecp%atoms(3)%local%nprim, mol_ecp%atoms(1)%local%nprim, &
                 "the second gold should have the same channels as the first")
      if (allocated(error)) return
      call check(error, close_enough(mol_ecp%atoms(3)%local%exponents(1), &
                                     mol_ecp%atoms(1)%local%exponents(1)), &
                 "the second gold's exponents should match the first's")

      call mol_ecp%destroy()
   end subroutine test_molecular

   subroutine test_core_sum(error)
      !! Core electrons sum over atoms, not over elements
      !!
      !! Two golds replace 120 electrons, not 60. Summing per element would
      !! leave the SCF with 60 electrons too many and still converge.
      type(error_type), allocatable, intent(out) :: error
      type(molecular_ecp_type) :: mol_ecp
      type(error_t) :: read_error
      character(len=2) :: symbols(3)

      symbols = ["Au", "Au", "O "]
      call build_molecular_ecp_json(ECP_FILE, symbols, mol_ecp, read_error)
      call check(error,.not. read_error%has_error(), read_error%get_message())
      if (allocated(error)) return

      call check(error, mol_ecp%core_electrons(), 120, &
                 "two gold ECPs should replace 120 electrons between them")
      if (allocated(error)) return
      call check(error, mol_ecp%any_ecp(), "the molecule should report carrying ECPs")

      call mol_ecp%destroy()
   end subroutine test_core_sum

   subroutine test_missing_file(error)
      type(error_type), allocatable, intent(out) :: error
      type(atomic_ecp_type) :: ecp
      type(error_t) :: read_error

      call read_json_ecp_element("../basis_sets/no-such-ecp.json", "Au", ecp, read_error)
      call check(error, read_error%has_error(), &
                 "a missing ECP file should be an error")
   end subroutine test_missing_file

   ! ---- helpers -------------------------------------------------------------

   pure function projected_of(ecp, ang_mom) result(index)
      !! Which projected channel carries a given l, or 0
      type(atomic_ecp_type), intent(in) :: ecp
      integer, intent(in) :: ang_mom
      integer :: index
      integer :: i

      index = 0
      do i = 1, ecp%n_projected
         if (ecp%projected(i)%ang_mom == ang_mom) then
            index = i
            return
         end if
      end do
   end function projected_of

   pure function close_enough(a, b) result(same)
      real(dp), intent(in) :: a, b
      logical :: same
      same = abs(a - b) <= 1.0e-10_dp*max(1.0_dp, abs(b))
   end function close_enough

end module test_mqc_ecp_reader

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_ecp_reader, only: collect_mqc_ecp_reader_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_ecp_reader", collect_mqc_ecp_reader_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
