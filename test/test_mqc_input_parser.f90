module test_mqc_input_parser
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_input_parser, only: input_config_t, read_input_file
   implicit none
   private
   public :: collect_mqc_input_parser_tests

contains

   !> Collect all exported unit tests
   subroutine collect_mqc_input_parser_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("read_valid_input", test_read_valid_input), &
                  new_unittest("read_input_with_comments", test_read_with_comments), &
                  new_unittest("read_input_with_quotes", test_read_with_quotes), &
                  new_unittest("read_input_missing_file", test_missing_file), &
                  new_unittest("read_input_missing_geom", test_missing_geom), &
                  new_unittest("read_input_missing_monomer", test_missing_monomer), &
                  new_unittest("read_input_custom_nlevel", test_custom_nlevel), &
                  new_unittest("read_input_invalid_nlevel", test_invalid_nlevel), &
                  new_unittest("read_input_with_gfn1", test_method_gfn1), &
                  new_unittest("read_input_with_gfn2", test_method_gfn2), &
                  new_unittest("read_input_invalid_method", test_invalid_method), &
                  new_unittest("read_input_default_method", test_default_method), &
                  new_unittest("config_destroy", test_config_destroy) &
                  ]
   end subroutine collect_mqc_input_parser_tests

   subroutine test_read_valid_input(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      ! Create test input file
      call create_test_input_basic()

      ! Read the input file
      call read_input_file("test_input_basic.inp", config, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to read input: "//errmsg)
         call cleanup_test_files()
         return
      end if

      call check(error, allocated(config%geom_file), "geom_file should be allocated")
      if (allocated(error)) then
         call config%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%geom_file), "water_cluster.xyz", &
                 "geom_file should be water_cluster.xyz")
      if (allocated(error)) then
         call config%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%monomer_file), "water.xyz", &
                 "monomer_file should be water.xyz")
      if (allocated(error)) then
         call config%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, config%nlevel, 1, "nlevel should default to 1")
      if (allocated(error)) then
         call config%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%method), "gfn2", "method should default to gfn2")

      call config%destroy()
      call cleanup_test_files()
   end subroutine test_read_valid_input

   subroutine test_read_with_comments(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_with_comments()

      call read_input_file("test_input_comments.inp", config, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to read input with comments: "//errmsg)
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%geom_file), "system.xyz", &
                 "Should read geom_file despite comments")

      call config%destroy()
      call cleanup_test_files()
   end subroutine test_read_with_comments

   subroutine test_read_with_quotes(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_with_quotes()

      call read_input_file("test_input_quotes.inp", config, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to read input with quotes: "//errmsg)
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%geom_file), "path/to/geom.xyz", &
                 "Should remove quotes from geom_file")
      if (allocated(error)) then
         call config%destroy()
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%monomer_file), "path/to/monomer.xyz", &
                 "Should remove quotes from monomer_file")

      call config%destroy()
      call cleanup_test_files()
   end subroutine test_read_with_quotes

   subroutine test_missing_file(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      ! Try to read non-existent file
      call read_input_file("nonexistent.inp", config, stat, errmsg)

      call check(error, stat /= 0, "Should fail with non-existent file")
   end subroutine test_missing_file

   subroutine test_missing_geom(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_missing_geom()

      call read_input_file("test_input_missing_geom.inp", config, stat, errmsg)

      call check(error, stat /= 0, "Should fail when geom field is missing")

      call cleanup_test_files()
   end subroutine test_missing_geom

   subroutine test_missing_monomer(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_missing_monomer()

      call read_input_file("test_input_missing_monomer.inp", config, stat, errmsg)

      call check(error, stat /= 0, "Should fail when monomer_symbols field is missing")

      call cleanup_test_files()
   end subroutine test_missing_monomer

   subroutine test_custom_nlevel(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_with_nlevel()

      call read_input_file("test_input_nlevel.inp", config, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to read input with nlevel: "//errmsg)
         call cleanup_test_files()
         return
      end if

      call check(error, config%nlevel, 3, "nlevel should be 3")

      call config%destroy()
      call cleanup_test_files()
   end subroutine test_custom_nlevel

   subroutine test_invalid_nlevel(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_invalid_nlevel()

      call read_input_file("test_input_invalid_nlevel.inp", config, stat, errmsg)

      call check(error, stat /= 0, "Should fail with negative nlevel")

      call cleanup_test_files()
   end subroutine test_invalid_nlevel

   subroutine test_config_destroy(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config

      config%geom_file = "test.xyz"
      config%monomer_file = "monomer.xyz"
      config%method = "gfn1"
      config%nlevel = 5

      call config%destroy()

      call check(error,.not. allocated(config%geom_file), &
                 "geom_file should be deallocated")
      if (allocated(error)) return

      call check(error,.not. allocated(config%monomer_file), &
                 "monomer_file should be deallocated")
      if (allocated(error)) return

      call check(error,.not. allocated(config%method), &
                 "method should be deallocated")
   end subroutine test_config_destroy

   subroutine test_method_gfn1(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_with_gfn1()

      call read_input_file("test_input_gfn1.inp", config, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to read input with gfn1: "//errmsg)
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%method), "gfn1", "method should be gfn1")

      call config%destroy()
      call cleanup_test_files()
   end subroutine test_method_gfn1

   subroutine test_method_gfn2(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_with_gfn2()

      call read_input_file("test_input_gfn2.inp", config, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to read input with gfn2: "//errmsg)
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%method), "gfn2", "method should be gfn2")

      call config%destroy()
      call cleanup_test_files()
   end subroutine test_method_gfn2

   subroutine test_invalid_method(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_invalid_method()

      call read_input_file("test_input_invalid_method.inp", config, stat, errmsg)

      call check(error, stat /= 0, "Should fail with invalid method")

      call cleanup_test_files()
   end subroutine test_invalid_method

   subroutine test_default_method(error)
      type(error_type), allocatable, intent(out) :: error
      type(input_config_t) :: config
      integer :: stat
      character(len=:), allocatable :: errmsg

      call create_test_input_basic()  ! No method specified

      call read_input_file("test_input_basic.inp", config, stat, errmsg)

      if (stat /= 0) then
         call check(error, .false., "Failed to read basic input: "//errmsg)
         call cleanup_test_files()
         return
      end if

      call check(error, trim(config%method), "gfn2", "method should default to gfn2")

      call config%destroy()
      call cleanup_test_files()
   end subroutine test_default_method

   ! Helper subroutines to create test files

   subroutine create_test_input_basic()
      integer :: unit
      open (newunit=unit, file="test_input_basic.inp", status="replace")
      write (unit, '(a)') "geom=water_cluster.xyz"
      write (unit, '(a)') "monomer_symbols=water.xyz"
      close (unit)
   end subroutine create_test_input_basic

   subroutine create_test_input_with_comments()
      integer :: unit
      open (newunit=unit, file="test_input_comments.inp", status="replace")
      write (unit, '(a)') "# This is a comment"
      write (unit, '(a)') "geom=system.xyz"
      write (unit, '(a)') "! Another comment style"
      write (unit, '(a)') "monomer_symbols=mono.xyz"
      write (unit, '(a)') ""
      write (unit, '(a)') "# nlevel=999  # this should be ignored"
      close (unit)
   end subroutine create_test_input_with_comments

   subroutine create_test_input_with_quotes()
      integer :: unit
      open (newunit=unit, file="test_input_quotes.inp", status="replace")
      write (unit, '(a)') 'geom="path/to/geom.xyz"'
      write (unit, '(a)') "monomer_symbols='path/to/monomer.xyz'"
      close (unit)
   end subroutine create_test_input_with_quotes

   subroutine create_test_input_missing_geom()
      integer :: unit
      open (newunit=unit, file="test_input_missing_geom.inp", status="replace")
      write (unit, '(a)') "monomer_symbols=water.xyz"
      close (unit)
   end subroutine create_test_input_missing_geom

   subroutine create_test_input_missing_monomer()
      integer :: unit
      open (newunit=unit, file="test_input_missing_monomer.inp", status="replace")
      write (unit, '(a)') "geom=water_cluster.xyz"
      close (unit)
   end subroutine create_test_input_missing_monomer

   subroutine create_test_input_with_nlevel()
      integer :: unit
      open (newunit=unit, file="test_input_nlevel.inp", status="replace")
      write (unit, '(a)') "geom=cluster.xyz"
      write (unit, '(a)') "monomer_symbols=monomer.xyz"
      write (unit, '(a)') "nlevel=3"
      close (unit)
   end subroutine create_test_input_with_nlevel

   subroutine create_test_input_invalid_nlevel()
      integer :: unit
      open (newunit=unit, file="test_input_invalid_nlevel.inp", status="replace")
      write (unit, '(a)') "geom=cluster.xyz"
      write (unit, '(a)') "monomer_symbols=monomer.xyz"
      write (unit, '(a)') "nlevel=-1"
      close (unit)
   end subroutine create_test_input_invalid_nlevel

   subroutine create_test_input_with_gfn1()
      integer :: unit
      open (newunit=unit, file="test_input_gfn1.inp", status="replace")
      write (unit, '(a)') "geom=cluster.xyz"
      write (unit, '(a)') "monomer_symbols=monomer.xyz"
      write (unit, '(a)') "method=gfn1"
      close (unit)
   end subroutine create_test_input_with_gfn1

   subroutine create_test_input_with_gfn2()
      integer :: unit
      open (newunit=unit, file="test_input_gfn2.inp", status="replace")
      write (unit, '(a)') "geom=cluster.xyz"
      write (unit, '(a)') "monomer_symbols=monomer.xyz"
      write (unit, '(a)') "method=gfn2"
      close (unit)
   end subroutine create_test_input_with_gfn2

   subroutine create_test_input_invalid_method()
      integer :: unit
      open (newunit=unit, file="test_input_invalid_method.inp", status="replace")
      write (unit, '(a)') "geom=cluster.xyz"
      write (unit, '(a)') "monomer_symbols=monomer.xyz"
      write (unit, '(a)') "method=gfn0"
      close (unit)
   end subroutine create_test_input_invalid_method

   subroutine cleanup_test_files()
      call delete_file("test_input_basic.inp")
      call delete_file("test_input_comments.inp")
      call delete_file("test_input_quotes.inp")
      call delete_file("test_input_missing_geom.inp")
      call delete_file("test_input_missing_monomer.inp")
      call delete_file("test_input_nlevel.inp")
      call delete_file("test_input_invalid_nlevel.inp")
      call delete_file("test_input_gfn1.inp")
      call delete_file("test_input_gfn2.inp")
      call delete_file("test_input_invalid_method.inp")
   end subroutine cleanup_test_files

   subroutine delete_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit, stat
      open (newunit=unit, file=filename, status="old", iostat=stat)
      if (stat == 0) close (unit, status="delete")
   end subroutine delete_file

end module test_mqc_input_parser

program tester
   use iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_input_parser, only: collect_mqc_input_parser_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_input_parser", collect_mqc_input_parser_tests) &
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
