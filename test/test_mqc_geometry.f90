module test_mqc_geometry
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_geometry, only: geometry_type
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_geometry_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_geometry_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("geometry_initialization", test_geometry_init), &
                  new_unittest("geometry_destroy", test_geometry_destroy) &
                  ]
   end subroutine collect_mqc_geometry_tests

   subroutine test_geometry_init(error)
      type(error_type), allocatable, intent(out) :: error
      type(geometry_type) :: geom

      ! Initialize a simple geometry
      geom%natoms = 2
      allocate (character(len=2) :: geom%elements(2))
      allocate (geom%coords(3, 2))
      allocate (character(len=20) :: geom%comment)

      geom%elements(1) = "O"
      geom%elements(2) = "H"
      geom%coords(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]
      geom%coords(:, 2) = [1.0_dp, 0.0_dp, 0.0_dp]
      geom%comment = "Test molecule"

      ! Check that everything is allocated and set correctly
      call check(error, geom%natoms, 2, "natoms should be 2")
      if (allocated(error)) return

      call check(error, allocated(geom%elements), "elements should be allocated")
      if (allocated(error)) return

      call check(error, allocated(geom%coords), "coords should be allocated")
      if (allocated(error)) return

      call check(error, allocated(geom%comment), "comment should be allocated")
      if (allocated(error)) return

      call check(error, trim(geom%elements(1)), "O", "first element should be O")
      if (allocated(error)) return

      call check(error, abs(geom%coords(1, 2) - 1.0_dp) < 1.0e-10_dp, &
                 "second atom x-coordinate should be 1.0")

      call geom%destroy()
   end subroutine test_geometry_init

   subroutine test_geometry_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(geometry_type) :: geom

      ! Create and populate geometry
      geom%natoms = 3
      allocate (character(len=2) :: geom%elements(3))
      allocate (geom%coords(3, 3))
      allocate (character(len=30) :: geom%comment)

      ! Destroy it
      call geom%destroy()

      ! Check that everything is deallocated
      call check(error,.not. allocated(geom%elements), &
                 "elements should be deallocated")
      if (allocated(error)) return

      call check(error,.not. allocated(geom%coords), &
                 "coords should be deallocated")
      if (allocated(error)) return

      call check(error,.not. allocated(geom%comment), &
                 "comment should be deallocated")
      if (allocated(error)) return

      call check(error, geom%natoms == 0, "natoms should be reset to 0")
   end subroutine test_geometry_destroy

end module test_mqc_geometry

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_geometry, only: collect_mqc_geometry_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_geometry", collect_mqc_geometry_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if

end program tester
