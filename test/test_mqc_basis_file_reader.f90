module test_mqc_basis_file_reader
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_basis_file_reader, only: basis_file_t, open_basis_file, extract_element
   use mqc_basis_reader, only: parse_element_basis
   use mqc_cgto, only: atomic_basis_type
   use mqc_error, only: error_t
   implicit none
   private
   public :: collect_mqc_basis_file_reader_tests

contains

   !> Find the basis_sets directory, trying multiple possible paths
   function find_basis_file() result(path)
      character(len=:), allocatable :: path
      logical :: exists
      character(len=*), parameter :: basis_file = "6-31G.txt"

      ! Try different possible paths (all same length for array constructor)
      character(len=50), parameter :: paths(4) = [ &
                                      character(len=50) :: &
                                      "../basis_sets/"//basis_file, &  ! CMake working directory
                                      "basis_sets/"//basis_file, &  ! fpm from root
                                      "../../basis_sets/"//basis_file, &  ! fpm from build dir
                                      "../../../basis_sets/"//basis_file]  ! deeply nested

      integer :: i

      do i = 1, size(paths)
         inquire (file=trim(paths(i)), exist=exists)
         if (exists) then
            path = trim(paths(i))
            return
         end if
      end do

      ! Not found, return the default
      path = "../basis_sets/"//basis_file
   end function find_basis_file

   !> Collect all exported unit tests
   subroutine collect_mqc_basis_file_reader_tests(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("open_6-31G_file", test_open_basis_file), &
                  new_unittest("extract_hydrogen", test_extract_hydrogen), &
                  new_unittest("extract_carbon", test_extract_carbon), &
                  new_unittest("parse_extracted_hydrogen", test_parse_extracted_hydrogen) &
                  ]
   end subroutine collect_mqc_basis_file_reader_tests

   subroutine test_open_basis_file(error)
      type(error_type), allocatable, intent(out) :: error
      type(basis_file_t) :: basis_file
      character(len=:), allocatable :: path_to_basis
      logical :: file_exists
      type(error_t) :: file_error

      ! Find the basis file
      path_to_basis = find_basis_file()
      inquire (file=path_to_basis, exist=file_exists)

      if (.not. file_exists) then
         call check(error, .false., "Basis file not found: "//path_to_basis)
         return
      end if

      call open_basis_file(basis_file, path_to_basis, file_error)
      call check(error,.not. file_error%has_error(), "Should open basis file without error")
      if (allocated(error)) return

      call check(error, len(basis_file%data_section) > 0, &
                 "Basis file should contain data")
      if (allocated(error)) return

      ! data_section is extracted *after* $DATA marker, so check for element names instead
      call check(error, index(basis_file%data_section, "HYDROGEN") > 0 .or. &
                 index(basis_file%data_section, "CARBON") > 0, &
                 "Basis file should contain element data")
   end subroutine test_open_basis_file

   subroutine test_extract_hydrogen(error)
      type(error_type), allocatable, intent(out) :: error
      type(basis_file_t) :: basis_file
      character(len=:), allocatable :: h_content, path_to_basis
      logical :: file_exists
      type(error_t) :: file_error

      path_to_basis = find_basis_file()
      inquire (file=path_to_basis, exist=file_exists)
      if (.not. file_exists) then
         call check(error, .false., "Basis file not found: "//path_to_basis)
         return
      end if

      call open_basis_file(basis_file, path_to_basis, file_error)
      if (file_error%has_error()) then
         call check(error, .false., "Failed to open basis file")
         return
      end if

      call extract_element(basis_file, "HYDROGEN", h_content, file_error)
      if (file_error%has_error()) then
         call check(error, .false., "Failed to extract hydrogen")
         return
      end if
      print *, h_content

      call check(error, len(h_content) > 0, "Should extract hydrogen basis content")
      if (allocated(error)) return

      call check(error, index(h_content, "HYDROGEN") > 0, &
                 "Extracted content should contain element name")
      if (allocated(error)) return

      call check(error, index(h_content, "S") > 0, &
                 "Hydrogen basis should contain S shell")
   end subroutine test_extract_hydrogen

   subroutine test_extract_carbon(error)
      type(error_type), allocatable, intent(out) :: error
      type(basis_file_t) :: basis_file
      character(len=:), allocatable :: c_content, path_to_basis
      logical :: file_exists
      type(error_t) :: file_error

      path_to_basis = find_basis_file()
      inquire (file=path_to_basis, exist=file_exists)
      if (.not. file_exists) then
         call check(error, .false., "Basis file not found: "//path_to_basis)
         return
      end if

      call open_basis_file(basis_file, path_to_basis, file_error)
      if (file_error%has_error()) then
         call check(error, .false., "Failed to open basis file")
         return
      end if

      call extract_element(basis_file, "CARBON", c_content, file_error)
      if (file_error%has_error()) then
         call check(error, .false., "Failed to extract carbon")
         return
      end if

      call check(error, len(c_content) > 0, "Should extract carbon basis content")
      if (allocated(error)) return

      call check(error, index(c_content, "CARBON") > 0, &
                 "Extracted content should contain element name")
      if (allocated(error)) return

      ! Carbon should have both S and L shells in 6-31G
      call check(error, index(c_content, "S") > 0, &
                 "Carbon basis should contain S shell")
      if (allocated(error)) return

      call check(error, index(c_content, "L") > 0, &
                 "Carbon basis should contain L shell")
   end subroutine test_extract_carbon

   subroutine test_parse_extracted_hydrogen(error)
      type(error_type), allocatable, intent(out) :: error
      type(basis_file_t) :: basis_file
      type(atomic_basis_type) :: h_basis
      character(len=:), allocatable :: h_content, path_to_basis
      type(error_t) :: parse_error, file_error
      logical :: file_exists

      path_to_basis = find_basis_file()
      inquire (file=path_to_basis, exist=file_exists)
      if (.not. file_exists) then
         call check(error, .false., "Basis file not found: "//path_to_basis)
         return
      end if

      call open_basis_file(basis_file, path_to_basis, file_error)
      if (file_error%has_error()) then
         call check(error, .false., "Failed to open basis file")
         return
      end if

      call extract_element(basis_file, "HYDROGEN", h_content, file_error)
      if (file_error%has_error()) then
         call check(error, .false., "Failed to extract hydrogen")
         return
      end if

      ! Parse the extracted content
      call parse_element_basis(h_content, "HYDROGEN", h_basis, parse_error)

      call check(error,.not. parse_error%has_error(), "Should parse extracted hydrogen basis")
      if (allocated(error)) return

      call check(error, trim(h_basis%element), "HYDROGEN", &
                 "Parsed basis should be for hydrogen")
      if (allocated(error)) return

      call check(error, h_basis%nshells > 0, &
                 "Hydrogen should have at least one shell")
      if (allocated(error)) return

      call h_basis%destroy()
   end subroutine test_parse_extracted_hydrogen

end module test_mqc_basis_file_reader

program tester_mqc_basis_file_reader
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_basis_file_reader, only: collect_mqc_basis_file_reader_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_basis_file_reader", collect_mqc_basis_file_reader_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_basis_file_reader
