module test_mqc_capi_run
   !! Tests the guard rails `mqc_run` puts up before it ever needs MPI or a
   !! quantum-chemistry backend.
   !!
   !! A real run needs a live session and tblite, which is the MPI validation
   !! harness's job. But `mqc_run` refuses three ways before it gets there, and
   !! all three are reachable serially: a null system handle, a partition that
   !! cuts bonds nobody declared, and a call made with no session started. The
   !! second is the one that matters -- an undeclared cut is exactly the mistake
   !! that would otherwise fragment a molecule into uncapped radicals and return
   !! a confident wrong number -- so it is checked here rather than trusted.
   use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_int, c_double, &
                                                                             c_char, c_null_char
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_capi_run, only: mqc_run, mqc_run_last_error
   use mqc_capi_system, only: mqc_system_new, mqc_system_free, &
                              mqc_system_set_geometry, mqc_system_set_monomers, &
                              mqc_system_set_bonds
   implicit none
   private
   public :: collect_mqc_capi_run_tests

   integer(c_int), parameter :: MQC_FAIL = 1
   integer(c_int), parameter :: MQC_BAD_HANDLE = 2

contains

   subroutine collect_mqc_capi_run_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("null_system_is_refused", test_null_system), &
                  new_unittest("undeclared_cut_is_refused", test_undeclared_cut), &
                  new_unittest("no_session_is_refused", test_no_session) &
                  ]
   end subroutine collect_mqc_capi_run_tests

   subroutine test_null_system(error)
      !! A null system handle is caught first of all
      type(error_type), allocatable, intent(out) :: error
      integer(c_int) :: status
      real(c_double) :: energy
      character(kind=c_char), allocatable :: settings(:), label(:)

      settings = to_cchars('{}')
      label = to_cchars('t')
      status = mqc_run(c_null_ptr, c_null_ptr, int(size(settings), c_int), settings, &
                       int(size(label), c_int), label, 0_c_int, energy)
      call check(error, int(status), int(MQC_BAD_HANDLE), "a null system is a bad handle")
      if (allocated(error)) return
      call check(error, index(last_run_error(), "null system") > 0, &
                 "the message should name the null system, got: "//last_run_error())
   end subroutine test_null_system

   subroutine test_undeclared_cut(error)
      !! A partition that severs an undeclared bond is refused before any run
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      real(c_double) :: energy
      character(kind=c_char), allocatable :: settings(:), label(:)
      integer(c_int) :: znum(2), sizes(2), atoms(2), charge(2), mult(2)
      real(c_double) :: coord(6)

      ! Two oxygens 1.2 A apart, split one per monomer, with no bonds declared.
      handle = mqc_system_new()
      znum = [8, 8]
      coord = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
               1.2_c_double, 0.0_c_double, 0.0_c_double]
      status = mqc_system_set_geometry(handle, 2_c_int, znum, coord, 0_c_int, 1_c_int)
      sizes = [1, 1]; atoms = [0, 1]; charge = [0, 0]; mult = [1, 1]
      status = mqc_system_set_monomers(handle, 2_c_int, 1_c_int, sizes, atoms, charge, mult)

      settings = to_cchars('{}')
      label = to_cchars('t')
      status = mqc_run(handle, c_null_ptr, int(size(settings), c_int), settings, &
                       int(size(label), c_int), label, 0_c_int, energy)
      call check(error, int(status), int(MQC_FAIL), "an undeclared cut is refused")
      if (allocated(error)) return
      call check(error, index(last_run_error(), "never declared") > 0, &
                 "the message should explain the undeclared cut, got: "//last_run_error())

      call mqc_system_free(handle)
   end subroutine test_undeclared_cut

   subroutine test_no_session(error)
      !! A fully declared system still cannot run without a session started
      type(error_type), allocatable, intent(out) :: error
      type(c_ptr) :: handle
      integer(c_int) :: status
      real(c_double) :: energy
      character(kind=c_char), allocatable :: settings(:), label(:)
      integer(c_int) :: znum(3), sizes(1), atoms(3), charge(1), mult(1), dummy(1)

      ! A single water, one monomer, bonds declared (none to cut): valid input.
      handle = mqc_system_new()
      znum = [8, 1, 1]
      status = mqc_system_set_geometry(handle, 3_c_int, znum, water_xyz(), 0_c_int, 1_c_int)
      sizes = [3]; atoms = [0, 1, 2]; charge = [0]; mult = [1]
      status = mqc_system_set_monomers(handle, 1_c_int, 3_c_int, sizes, atoms, charge, mult)
      dummy = 0
      status = mqc_system_set_bonds(handle, 0_c_int, dummy, dummy, dummy, dummy)

      settings = to_cchars('{}')
      label = to_cchars('t')
      status = mqc_run(handle, c_null_ptr, int(size(settings), c_int), settings, &
                       int(size(label), c_int), label, 0_c_int, energy)
      call check(error, int(status), int(MQC_FAIL), "no session means no run")
      if (allocated(error)) return
      call check(error, index(last_run_error(), "session") > 0, &
                 "the message should point at the missing session, got: "//last_run_error())

      call mqc_system_free(handle)
   end subroutine test_no_session

   ! ---- helpers ------------------------------------------------------------

   pure function water_xyz() result(xyz)
      !! A single water, Angstrom, x,y,z contiguous per atom
      real(c_double) :: xyz(9)
      xyz = [0.0_c_double, 0.0_c_double, 0.0_c_double, &
             0.757_c_double, 0.586_c_double, 0.0_c_double, &
             -0.757_c_double, 0.586_c_double, 0.0_c_double]
   end function water_xyz

   function to_cchars(s) result(arr)
      !! A Fortran string as the C character array the entry points take
      character(len=*), intent(in) :: s
      character(kind=c_char), allocatable :: arr(:)
      integer :: i

      allocate (arr(len(s)))
      do i = 1, len(s)
         arr(i) = s(i:i)
      end do
   end function to_cchars

   function last_run_error() result(message)
      !! Pull the most recent run error out through the C string boundary
      character(len=:), allocatable :: message
      integer, parameter :: BUFLEN = 512
      character(kind=c_char) :: buffer(BUFLEN)
      integer :: i

      buffer = c_null_char
      call mqc_run_last_error(int(BUFLEN, c_int), buffer)
      message = ""
      do i = 1, BUFLEN
         if (buffer(i) == c_null_char) exit
         message = message//buffer(i)
      end do
   end function last_run_error

end module test_mqc_capi_run

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_capi_run, only: collect_mqc_capi_run_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_capi_run", collect_mqc_capi_run_tests) &
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
