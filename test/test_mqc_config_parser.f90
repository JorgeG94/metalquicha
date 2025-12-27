!! Unit tests for mqc_config_parser module
program test_mqc_config_parser
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mqc_config_parser, only: mqc_config_t, read_mqc_file
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
   use mqc_calc_types, only: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT
   use pic_test_helpers, only: is_equal
   implicit none
   integer :: stat
   type(unittest_type), allocatable :: testsuite(:)
   character(len=:), allocatable :: errmsg

   ! Build test suite
   testsuite = [ &
               new_unittest("parse_minimal", test_parse_minimal), &
               new_unittest("parse_with_fragments", test_parse_with_fragments), &
               new_unittest("parse_with_connectivity", test_parse_with_connectivity), &
               new_unittest("parse_no_fragments", test_parse_no_fragments), &
               new_unittest("parse_method_xtb", test_parse_method_xtb), &
               new_unittest("parse_system_log_level", test_parse_system_log_level), &
               new_unittest("error_missing_schema", test_error_missing_schema), &
               new_unittest("error_missing_geometry", test_error_missing_geometry) &
               ]

   call run_testsuite(testsuite, stat, errmsg)

   if (stat /= 0) then
      write (*, '(A)') "Test suite failed: "//errmsg
      stop 1
   else
      write (*, '(A)') "All tests passed!"
   end if

contains

   subroutine test_parse_minimal(error)
      !! Test parsing minimal .mqc file with required sections only
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      integer :: stat, unit
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_file = "test_minimal.mqc"

      ! Create minimal test file
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%model'
      write (unit, '(A)') 'method = XTB-GFN2'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%driver'
      write (unit, '(A)') 'type = Energy'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%structure'
      write (unit, '(A)') 'charge = 0'
      write (unit, '(A)') 'multiplicity = 1'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '2'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'H 0.7 0.0 0.0'
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, stat, errmsg)

      call check(error, stat, 0, "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%schema_name, "mqc-frag", "schema_name should be 'mqc-frag'")
      if (allocated(error)) return

      call check(error, config%schema_version, "1.0", "schema_version should be '1.0'")
      if (allocated(error)) return

      call check(error, config%index_base, 0, "index_base should be 0")
      if (allocated(error)) return

      call check(error, config%units, "angstrom", "units should be 'angstrom'")
      if (allocated(error)) return

      call check(error, config%method, METHOD_TYPE_GFN2, "method should be GFN2")
      if (allocated(error)) return

      call check(error, config%calc_type, CALC_TYPE_ENERGY, "calc_type should be ENERGY")
      if (allocated(error)) return

      call check(error, config%charge, 0, "charge should be 0")
      if (allocated(error)) return

      call check(error, config%multiplicity, 1, "multiplicity should be 1")
      if (allocated(error)) return

      call check(error, config%geometry%natoms, 2, "natoms should be 2")
      if (allocated(error)) return

      call check(error, config%geometry%elements(1), "H", "First element should be H")
      if (allocated(error)) return

      call check(error, is_equal(config%geometry%coords(1, 1), 0.0d0), &
                 "First atom x-coordinate should be 0.0")
      if (allocated(error)) return

      call check(error, is_equal(config%geometry%coords(1, 2), 0.7d0), &
                 "Second atom x-coordinate should be 0.7")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_minimal

   subroutine test_parse_with_fragments(error)
      !! Test parsing .mqc file with fragments
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      integer :: stat, unit
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_file = "test_fragments.mqc"

      ! Create test file with fragments
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%model'
      write (unit, '(A)') 'method = XTB-GFN1'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%driver'
      write (unit, '(A)') 'type = Gradient'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%structure'
      write (unit, '(A)') 'charge = 1'
      write (unit, '(A)') 'multiplicity = 2'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '4'
      write (unit, '(A)') 'Water dimer'
      write (unit, '(A)') 'O 0.0 0.0 0.0'
      write (unit, '(A)') 'H 1.0 0.0 0.0'
      write (unit, '(A)') 'O 3.0 0.0 0.0'
      write (unit, '(A)') 'H 4.0 0.0 0.0'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragments'
      write (unit, '(A)') 'nfrag = 2'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragment'
      write (unit, '(A)') 'charge = 0'
      write (unit, '(A)') 'multiplicity = 1'
      write (unit, '(A)') '%indices'
      write (unit, '(A)') '0 1'
      write (unit, '(A)') 'end'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragment'
      write (unit, '(A)') 'charge = 1'
      write (unit, '(A)') 'multiplicity = 2'
      write (unit, '(A)') '%indices'
      write (unit, '(A)') '2 3'
      write (unit, '(A)') 'end'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, stat, errmsg)

      call check(error, stat, 0, "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%method, METHOD_TYPE_GFN1, "method should be GFN1")
      if (allocated(error)) return

      call check(error, config%calc_type, CALC_TYPE_GRADIENT, "calc_type should be GRADIENT")
      if (allocated(error)) return

      call check(error, config%charge, 1, "charge should be 1")
      if (allocated(error)) return

      call check(error, config%multiplicity, 2, "multiplicity should be 2")
      if (allocated(error)) return

      call check(error, config%nfrag, 2, "nfrag should be 2")
      if (allocated(error)) return

      call check(error, config%fragments(1)%charge, 0, "Fragment 1 charge should be 0")
      if (allocated(error)) return

      call check(error, config%fragments(1)%multiplicity, 1, "Fragment 1 multiplicity should be 1")
      if (allocated(error)) return

      call check(error, size(config%fragments(1)%indices), 2, "Fragment 1 should have 2 indices")
      if (allocated(error)) return

      call check(error, config%fragments(1)%indices(1), 0, "Fragment 1 first index should be 0")
      if (allocated(error)) return

      call check(error, config%fragments(2)%charge, 1, "Fragment 2 charge should be 1")
      if (allocated(error)) return

      call check(error, config%fragments(2)%multiplicity, 2, "Fragment 2 multiplicity should be 2")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_with_fragments

   subroutine test_parse_with_connectivity(error)
      !! Test parsing .mqc file with connectivity
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      integer :: stat, unit
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_file = "test_connectivity.mqc"

      ! Create test file with connectivity
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '2'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'H 0.7 0.0 0.0'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%connectivity'
      write (unit, '(A)') 'nbonds = 1'
      write (unit, '(A)') ''
      write (unit, '(A)') '0 1 1 broken'
      write (unit, '(A)') ''
      write (unit, '(A)') 'nbroken = 1'
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, stat, errmsg)

      call check(error, stat, 0, "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%nbonds, 1, "nbonds should be 1")
      if (allocated(error)) return

      call check(error, config%nbroken, 1, "nbroken should be 1")
      if (allocated(error)) return

      call check(error, config%bonds(1)%atom_i, 0, "Bond atom_i should be 0")
      if (allocated(error)) return

      call check(error, config%bonds(1)%atom_j, 1, "Bond atom_j should be 1")
      if (allocated(error)) return

      call check(error, config%bonds(1)%order, 1, "Bond order should be 1")
      if (allocated(error)) return

      call check(error, config%bonds(1)%is_broken, .true., "Bond should be broken")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_with_connectivity

   subroutine test_parse_no_fragments(error)
      !! Test parsing .mqc file without fragments (unfragmented calculation)
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      integer :: stat, unit
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_file = "test_no_fragments.mqc"

      ! Create test file without fragments
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '1'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, stat, errmsg)

      call check(error, stat, 0, "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%nfrag, 0, "nfrag should be 0 (no fragments)")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_no_fragments

   subroutine test_parse_method_xtb(error)
      !! Test parsing XTB method strings
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      integer :: stat, unit
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_file = "test_method.mqc"

      ! Create test file
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%model'
      write (unit, '(A)') 'method = XTB-GFN1'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '1'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, stat, errmsg)

      call check(error, stat, 0, "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%method, METHOD_TYPE_GFN1, "method should be GFN1 from XTB-GFN1")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_method_xtb

   subroutine test_parse_system_log_level(error)
      !! Test parsing %system section with log_level
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      integer :: stat, unit
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_file = "test_system.mqc"

      ! Create test file with %system section
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%system'
      write (unit, '(A)') 'log_level = verbose'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '1'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, stat, errmsg)

      call check(error, stat, 0, "Parser should succeed")
      if (allocated(error)) return

      call check(error, allocated(config%log_level), "log_level should be allocated")
      if (allocated(error)) return

      call check(error, config%log_level, "verbose", "log_level should be 'verbose'")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_system_log_level

   subroutine test_error_missing_schema(error)
      !! Test error handling when schema section is missing
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      integer :: stat, unit
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_file = "test_missing_schema.mqc"

      ! Create test file without schema
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '1'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, stat, errmsg)

      call check(error, stat /= 0, "Parser should fail for missing schema")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_error_missing_schema

   subroutine test_error_missing_geometry(error)
      !! Test error handling when geometry section is missing
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      integer :: stat, unit
      character(len=:), allocatable :: errmsg
      character(len=*), parameter :: test_file = "test_missing_geometry.mqc"

      ! Create test file without geometry
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, stat, errmsg)

      call check(error, stat /= 0, "Parser should fail for missing geometry")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_error_missing_geometry

   subroutine run_testsuite(tests, stat, errmsg)
      !! Simple test runner
      type(unittest_type), intent(in) :: tests(:)
      integer, intent(out) :: stat
      character(len=:), allocatable, intent(out) :: errmsg

      type(error_type), allocatable :: error
      integer :: i, passed, failed

      stat = 0
      passed = 0
      failed = 0

      do i = 1, size(tests)
         write (*, '(A,I0,A,I0,A)', advance='no') "Running test ", i, "/", size(tests), ": "
         write (*, '(A)', advance='no') trim(tests(i)%name)//"..."

         call tests(i)%test(error)

         if (allocated(error)) then
            write (*, '(A)') " FAILED"
            write (*, '(A)') "  Error: "//error%message
            failed = failed + 1
            stat = 1
            if (.not. allocated(errmsg)) then
               errmsg = "Test '"//trim(tests(i)%name)//"' failed: "//error%message
            end if
            deallocate (error)
         else
            write (*, '(A)') " PASSED"
            passed = passed + 1
         end if
      end do

      write (*, '(/,A,I0,A,I0,A)') "Tests passed: ", passed, "/", size(tests), ""
      if (failed > 0) then
         write (*, '(A,I0)') "Tests failed: ", failed
      end if

   end subroutine run_testsuite

end program test_mqc_config_parser
