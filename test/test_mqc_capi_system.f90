module test_mqc_capi_system
   !! Tests the C entry points that describe a system to be fragmented.
   !!
   !! These are `bind(C)` procedures, but they are ordinary Fortran module
   !! procedures too, so the suite calls them the way C would -- opaque handle,
   !! `c_int` status, 0-based atom indices, coordinates in Angstrom -- and
   !! checks the contract the header promises rather than the Fortran internals
   !! behind it. The value is in the guard rails: every "set the geometry first",
   !! every out-of-range index, and the one failure worth failing hard over -- a
   !! bond that crosses a monomer boundary without being marked broken.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_double, &
                                                                             c_char, c_null_char
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_capi_system, only: mqc_system_new, mqc_system_free, &
                              mqc_system_set_geometry, mqc_system_set_monomers, &
                              mqc_system_set_bonds, mqc_system_n_atoms, &
                              mqc_system_n_monomers, mqc_system_n_bonds, &
                              mqc_system_bonds_declared, mqc_system_perceive_bonds, &
                              mqc_system_count_missing_bonds, mqc_system_auto_monomers, &
                              mqc_system_last_error
   implicit none
   private
   public :: collect_mqc_capi_system_tests

   integer(c_int), parameter :: MQC_OK = 0
   integer(c_int), parameter :: MQC_FAIL = 1
   integer(c_int), parameter :: MQC_BAD_HANDLE = 2

contains

   subroutine collect_mqc_capi_system_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("empty_handle_reports_zeros", test_empty_counts), &
                  new_unittest("null_handle_is_refused", test_null_handle), &
                  new_unittest("set_geometry_happy_path", test_set_geometry), &
                  new_unittest("set_geometry_validates", test_geometry_validation), &
                  new_unittest("monomers_need_geometry_first", test_monomers_order), &
                  new_unittest("set_monomers_happy_path", test_set_monomers), &
                  new_unittest("set_monomers_validates", test_monomer_validation), &
                  new_unittest("bonds_need_monomers_first", test_bonds_order), &
                  new_unittest("zero_bonds_is_a_declaration", test_zero_bonds), &
                  new_unittest("set_bonds_records_them", test_set_bonds), &
                  new_unittest("uncut_boundary_bond_is_refused", test_uncut_boundary), &
                  new_unittest("set_bonds_validates", test_bond_validation), &
                  new_unittest("perceive_bonds_finds_intramolecular", test_perceive), &
                  new_unittest("count_missing_flags_undeclared_cut", test_count_missing), &
                  new_unittest("auto_monomers_splits_a_cluster", test_auto_cluster), &
                  new_unittest("auto_monomers_refuses_one_molecule", test_auto_single), &
                  new_unittest("last_error_is_retrievable", test_last_error) &
                  ]
   end subroutine collect_mqc_capi_system_tests

   subroutine test_empty_counts(error)
      !! A fresh handle has nothing in it, and says so with zeros not junk
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle

      handle = mqc_system_new()
      call check(error, int(mqc_system_n_atoms(handle)), 0, "a new system has no atoms")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_monomers(handle)), 0, "no monomers yet")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_bonds(handle)), 0, "no bonds yet")
      if (allocated(error)) return
      call check(error, int(mqc_system_bonds_declared(handle)), 0, &
                 "bonds have not been declared on a fresh handle")

      call mqc_system_free(handle)
   end subroutine test_empty_counts

   subroutine test_null_handle(error)
      !! The queries answer -1 for a bad handle, the setters refuse it
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: status
      integer(c_int) :: znum(1), charge(1), mult(1), sizes(1), atoms(1)
      real(c_double) :: coord(3)

      call check(error, int(mqc_system_n_atoms(c_null_ptr)), -1, "null -> -1 atoms")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_monomers(c_null_ptr)), -1, "null -> -1 monomers")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_bonds(c_null_ptr)), -1, "null -> -1 bonds")
      if (allocated(error)) return
      call check(error, int(mqc_system_bonds_declared(c_null_ptr)), 0, &
                 "a null handle has declared nothing")
      if (allocated(error)) return

      znum = 8; charge = 0; mult = 1; sizes = 1; atoms = 0; coord = 0.0_c_double
      status = mqc_system_set_geometry(c_null_ptr, 1_c_int, znum, coord, 0_c_int, 1_c_int)
      call check(error, int(status), int(MQC_BAD_HANDLE), "null handle -> bad handle")
      if (allocated(error)) return
      status = mqc_system_set_monomers(c_null_ptr, 1_c_int, 1_c_int, sizes, atoms, charge, mult)
      call check(error, int(status), int(MQC_BAD_HANDLE), "null handle -> bad handle")
   end subroutine test_null_handle

   subroutine test_set_geometry(error)
      !! Two waters go in and the atom count comes back
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status

      handle = mqc_system_new()
      status = mqc_system_set_geometry(handle, 6_c_int, two_waters_z(), two_waters_xyz(), &
                                       0_c_int, 1_c_int)
      call check(error, int(status), int(MQC_OK), "a valid geometry is accepted")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_atoms(handle)), 6, "six atoms in two waters")

      call mqc_system_free(handle)
   end subroutine test_set_geometry

   subroutine test_geometry_validation(error)
      !! No atoms and a zero multiplicity are both refused
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: znum(1)
      real(c_double) :: coord(3)

      handle = mqc_system_new()
      znum = 1; coord = 0.0_c_double
      status = mqc_system_set_geometry(handle, 0_c_int, znum, coord, 0_c_int, 1_c_int)
      call check(error, int(status), int(MQC_FAIL), "zero atoms is an error")
      if (allocated(error)) return

      status = mqc_system_set_geometry(handle, 1_c_int, znum, coord, 0_c_int, 0_c_int)
      call check(error, int(status), int(MQC_FAIL), "multiplicity below 1 is an error")

      call mqc_system_free(handle)
   end subroutine test_geometry_validation

   subroutine test_monomers_order(error)
      !! Monomers cannot be set before there is a geometry to partition
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: sizes(1), atoms(1), charge(1), mult(1)

      handle = mqc_system_new()
      sizes = 1; atoms = 0; charge = 0; mult = 1
      status = mqc_system_set_monomers(handle, 1_c_int, 1_c_int, sizes, atoms, charge, mult)
      call check(error, int(status), int(MQC_FAIL), "set the geometry first")

      call mqc_system_free(handle)
   end subroutine test_monomers_order

   subroutine test_set_monomers(error)
      !! Two waters, two monomers, column per monomer
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: sizes(2), atoms(6), charge(2), mult(2)

      handle = mqc_system_new()
      status = mqc_system_set_geometry(handle, 6_c_int, two_waters_z(), two_waters_xyz(), &
                                       0_c_int, 1_c_int)
      call check(error, int(status), int(MQC_OK), "geometry accepted")
      if (allocated(error)) return

      sizes = [3, 3]
      atoms = [0, 1, 2, 3, 4, 5]   ! column per monomer, 0-based
      charge = [0, 0]
      mult = [1, 1]
      status = mqc_system_set_monomers(handle, 2_c_int, 3_c_int, sizes, atoms, charge, mult)
      call check(error, int(status), int(MQC_OK), "a valid partition is accepted")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_monomers(handle)), 2, "two monomers")

      call mqc_system_free(handle)
   end subroutine test_set_monomers

   subroutine test_monomer_validation(error)
      !! A size out of range and an atom index out of range are both refused
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: sizes(2), atoms(6), charge(2), mult(2)

      handle = mqc_system_new()
      status = mqc_system_set_geometry(handle, 6_c_int, two_waters_z(), two_waters_xyz(), &
                                       0_c_int, 1_c_int)
      charge = [0, 0]; mult = [1, 1]; atoms = [0, 1, 2, 3, 4, 5]

      sizes = [4, 3]   ! 4 > max_size = 3
      status = mqc_system_set_monomers(handle, 2_c_int, 3_c_int, sizes, atoms, charge, mult)
      call check(error, int(status), int(MQC_FAIL), "a size above max_size is refused")
      if (allocated(error)) return

      sizes = [3, 3]
      atoms = [0, 1, 2, 3, 4, 9]   ! 9 >= total_atoms = 6
      status = mqc_system_set_monomers(handle, 2_c_int, 3_c_int, sizes, atoms, charge, mult)
      call check(error, int(status), int(MQC_FAIL), "an out-of-range atom index is refused")

      call mqc_system_free(handle)
   end subroutine test_monomer_validation

   subroutine test_bonds_order(error)
      !! Bonds cannot be judged before the partition that decides which are cut
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: ai(1), aj(1), ord(1), brk(1)

      handle = mqc_system_new()
      status = mqc_system_set_geometry(handle, 6_c_int, two_waters_z(), two_waters_xyz(), &
                                       0_c_int, 1_c_int)
      ai = 0; aj = 1; ord = 1; brk = 0
      status = mqc_system_set_bonds(handle, 1_c_int, ai, aj, ord, brk)
      call check(error, int(status), int(MQC_FAIL), "set the monomers first")

      call mqc_system_free(handle)
   end subroutine test_bonds_order

   subroutine test_zero_bonds(error)
      !! "No bonds" is a statement, and it flips bonds_declared
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: dummy(1)

      handle = two_water_partition()
      dummy = 0
      status = mqc_system_set_bonds(handle, 0_c_int, dummy, dummy, dummy, dummy)
      call check(error, int(status), int(MQC_OK), "declaring zero bonds is allowed")
      if (allocated(error)) return
      call check(error, int(mqc_system_bonds_declared(handle)), 1, &
                 "even zero bonds counts as having declared them")

      call mqc_system_free(handle)
   end subroutine test_zero_bonds

   subroutine test_set_bonds(error)
      !! An intramolecular bond and a broken cross-boundary bond are recorded
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: ai(2), aj(2), ord(2), brk(2)

      handle = two_water_partition()
      ai = [0, 0]        ! O-H inside monomer 1, then O(1)-O(2) across the boundary
      aj = [1, 3]
      ord = [1, 1]
      brk = [0, 1]       ! the boundary-crossing one is marked broken
      status = mqc_system_set_bonds(handle, 2_c_int, ai, aj, ord, brk)
      call check(error, int(status), int(MQC_OK), "valid bonds are accepted")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_bonds(handle)), 2, "two bonds recorded")
      if (allocated(error)) return
      call check(error, int(mqc_system_bonds_declared(handle)), 1, "bonds now declared")

      call mqc_system_free(handle)
   end subroutine test_set_bonds

   subroutine test_uncut_boundary(error)
      !! The one failure worth failing hard over
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: ai(1), aj(1), ord(1), brk(1)

      handle = two_water_partition()
      ai = 0; aj = 3     ! atoms in different monomers
      ord = 1; brk = 0   ! ... but not marked broken
      status = mqc_system_set_bonds(handle, 1_c_int, ai, aj, ord, brk)
      call check(error, int(status), int(MQC_FAIL), &
                 "a boundary-crossing bond left unbroken must be refused")
      if (allocated(error)) return
      call check(error, int(mqc_system_bonds_declared(handle)), 0, &
                 "a refused set_bonds declares nothing")

      call mqc_system_free(handle)
   end subroutine test_uncut_boundary

   subroutine test_bond_validation(error)
      !! Out-of-range, self-bonded, and zero-order bonds are all refused
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: ai(1), aj(1), ord(1), brk(1)

      handle = two_water_partition()

      ai = 0; aj = 9; ord = 1; brk = 1   ! atom 9 does not exist
      status = mqc_system_set_bonds(handle, 1_c_int, ai, aj, ord, brk)
      call check(error, int(status), int(MQC_FAIL), "an out-of-range atom is refused")
      if (allocated(error)) return

      ai = 2; aj = 2; ord = 1; brk = 0   ! an atom bonded to itself
      status = mqc_system_set_bonds(handle, 1_c_int, ai, aj, ord, brk)
      call check(error, int(status), int(MQC_FAIL), "a self-bond is refused")
      if (allocated(error)) return

      ai = 0; aj = 1; ord = 0; brk = 0   ! order below 1
      status = mqc_system_set_bonds(handle, 1_c_int, ai, aj, ord, brk)
      call check(error, int(status), int(MQC_FAIL), "a zero bond order is refused")

      call mqc_system_free(handle)
   end subroutine test_bond_validation

   subroutine test_perceive(error)
      !! Perception finds the four intramolecular O-H bonds and no cut
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status

      handle = two_water_partition()
      status = mqc_system_perceive_bonds(handle, -1.0_c_double)   ! <= 0 asks for the default
      call check(error, int(status), int(MQC_OK), "perception succeeds")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_bonds(handle)), 4, "two O-H bonds per water")
      if (allocated(error)) return
      call check(error, int(mqc_system_bonds_declared(handle)), 1, &
                 "perception declares the bonds it found")
      if (allocated(error)) return
      call check(error, int(mqc_system_count_missing_bonds(handle, -1.0_c_double)), 0, &
                 "perception leaves nothing missing")

      call mqc_system_free(handle)
   end subroutine test_perceive

   subroutine test_count_missing(error)
      !! A cut nobody declared is counted; a bad handle answers -1
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: znum(2), sizes(2), atoms(2), charge(2), mult(2)
      real(c_double) :: coord(6)

      call check(error, int(mqc_system_count_missing_bonds(c_null_ptr, -1.0_c_double)), -1, &
                 "a bad handle cannot answer, so -1")
      if (allocated(error)) return

      ! Two oxygens 1.2 A apart -- bonded by perception -- split one per monomer.
      handle = mqc_system_new()
      znum = [8, 8]
      coord = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
               1.2_c_double, 0.0_c_double, 0.0_c_double]
      status = mqc_system_set_geometry(handle, 2_c_int, znum, coord, 0_c_int, 1_c_int)
      sizes = [1, 1]; atoms = [0, 1]; charge = [0, 0]; mult = [1, 1]
      status = mqc_system_set_monomers(handle, 2_c_int, 1_c_int, sizes, atoms, charge, mult)

      call check(error, int(mqc_system_count_missing_bonds(handle, -1.0_c_double)), 1, &
                 "the partition cuts one bond nobody declared")

      call mqc_system_free(handle)
   end subroutine test_count_missing

   subroutine test_auto_cluster(error)
      !! Two separate waters become two monomers on their own
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status

      handle = mqc_system_new()
      status = mqc_system_set_geometry(handle, 6_c_int, two_waters_z(), two_waters_xyz(), &
                                       0_c_int, 1_c_int)
      status = mqc_system_auto_monomers(handle, -1.0_c_double)
      call check(error, int(status), int(MQC_OK), "a cluster partitions automatically")
      if (allocated(error)) return
      call check(error, int(mqc_system_n_monomers(handle)), 2, "one monomer per molecule")

      call mqc_system_free(handle)
   end subroutine test_auto_cluster

   subroutine test_auto_single(error)
      !! A single connected molecule has no automatic partition, and says so
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: znum(3)
      real(c_double) :: coord(9)

      handle = mqc_system_new()
      znum = [8, 1, 1]
      coord = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
               0.757_c_double, 0.586_c_double, 0.0_c_double, &
               -0.757_c_double, 0.586_c_double, 0.0_c_double]
      status = mqc_system_set_geometry(handle, 3_c_int, znum, coord, 0_c_int, 1_c_int)
      status = mqc_system_auto_monomers(handle, -1.0_c_double)
      call check(error, int(status), int(MQC_FAIL), &
                 "one molecule cannot be partitioned automatically")

      call mqc_system_free(handle)
   end subroutine test_auto_single

   subroutine test_last_error(error)
      !! After a refusal, the message can be copied out as a C string
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: znum(1)
      real(c_double) :: coord(3)
      character(len=:), allocatable :: message

      handle = mqc_system_new()
      znum = 1; coord = 0.0_c_double
      status = mqc_system_set_geometry(handle, 0_c_int, znum, coord, 0_c_int, 1_c_int)
      call check(error, int(status), int(MQC_FAIL), "provoke a failure")
      if (allocated(error)) return

      message = last_system_error()
      call check(error, len_trim(message) > 0, "the failure left a message behind")
      if (allocated(error)) return
      call check(error, index(message, "at least one atom") > 0, &
                 "the message names the actual problem, got: "//message)

      call mqc_system_free(handle)
   end subroutine test_last_error

   ! ---- fixtures and helpers ----------------------------------------------

   function two_water_partition() result(handle)
      !! A handle carrying two waters already split one molecule per monomer
      type(c_ptr) :: handle
      integer(c_int) :: status
      integer(c_int) :: sizes(2), atoms(6), charge(2), mult(2)

      handle = mqc_system_new()
      status = mqc_system_set_geometry(handle, 6_c_int, two_waters_z(), two_waters_xyz(), &
                                       0_c_int, 1_c_int)
      sizes = [3, 3]
      atoms = [0, 1, 2, 3, 4, 5]
      charge = [0, 0]
      mult = [1, 1]
      status = mqc_system_set_monomers(handle, 2_c_int, 3_c_int, sizes, atoms, charge, mult)
   end function two_water_partition

   pure function two_waters_z() result(z)
      !! Atomic numbers for O,H,H,O,H,H
      integer(c_int) :: z(6)
      z = [8, 1, 1, 8, 1, 1]
   end function two_waters_z

   pure function two_waters_xyz() result(xyz)
      !! Two waters 3 A apart, Angstrom, x,y,z contiguous per atom
      real(c_double) :: xyz(18)
      xyz = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
             0.757_c_double, 0.586_c_double, 0.0_c_double, &
             -0.757_c_double, 0.586_c_double, 0.0_c_double, &
             3.0_c_double, 0.0_c_double, 0.0_c_double, &
             3.757_c_double, 0.586_c_double, 0.0_c_double, &
             2.243_c_double, 0.586_c_double, 0.0_c_double]
   end function two_waters_xyz

   function last_system_error() result(message)
      !! Pull the most recent system error out through the C string boundary
      character(len=:), allocatable :: message
      integer, parameter :: BUFLEN = 512
      character(kind=c_char) :: buffer(BUFLEN)
      integer :: i

      buffer = c_null_char
      call mqc_system_last_error(int(BUFLEN, c_int), buffer)
      message = ""
      do i = 1, BUFLEN
         if (buffer(i) == c_null_char) exit
         message = message//buffer(i)
      end do
   end function last_system_error

end module test_mqc_capi_system

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_capi_system, only: collect_mqc_capi_system_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_capi_system", collect_mqc_capi_system_tests) &
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
