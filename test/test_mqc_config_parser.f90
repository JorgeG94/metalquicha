!! Unit tests for mqc_config_parser module
program test_mqc_config_parser
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use mqc_config_parser, only: mqc_config_t, read_mqc_file
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
   use mqc_calc_types, only: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, DEFAULT_TEMPERATURE, &
                                       DEFAULT_PRESSURE, DEFAULT_AIMD_DT, &
                                       DEFAULT_AIMD_NSTEPS, DEFAULT_AIMD_TEMPERATURE, &
                                       DEFAULT_AIMD_OUTPUT_FREQ
   use mqc_error, only: error_t
   use pic_types, only: dp
   use pic_test_helpers, only: is_equal
   implicit none
   type(error_t) :: error
   type(unittest_type), allocatable :: testsuite(:)

   ! Build test suite
   testsuite = [ &
               new_unittest("parse_minimal", test_parse_minimal), &
               new_unittest("parse_with_fragments", test_parse_with_fragments), &
               new_unittest("parse_with_connectivity", test_parse_with_connectivity), &
               new_unittest("parse_no_fragments", test_parse_no_fragments), &
               new_unittest("parse_method_xtb", test_parse_method_xtb), &
               new_unittest("parse_system_log_level", test_parse_system_log_level), &
               new_unittest("parse_with_comments", test_parse_with_comments), &
               new_unittest("parse_hessian_section", test_parse_hessian_section), &
               new_unittest("parse_hessian_defaults", test_parse_hessian_defaults), &
               new_unittest("parse_aimd_section", test_parse_aimd_section), &
               new_unittest("parse_fragmentation_section", test_parse_fragmentation_section), &
               new_unittest("parse_fragmentation_cutoffs", test_parse_fragmentation_cutoffs), &
               new_unittest("parse_fragmentation_global_groups", test_parse_fragmentation_global_groups), &
               new_unittest("parse_fragmentation_nodes_per_group", test_parse_fragmentation_nodes_per_group), &
               new_unittest("error_invalid_global_groups", test_error_invalid_global_groups), &
               new_unittest("parse_xtb_solvation", test_parse_xtb_solvation), &
               new_unittest("error_missing_schema", test_error_missing_schema), &
               new_unittest("error_missing_geometry", test_error_missing_geometry) &
               ]

   call run_testsuite(testsuite, error%code, error%message)

   if (error%has_error()) then
      write (*, '(A)') "Test suite failed: "//error%get_message()
      stop 1
   else
      write (*, '(A)') "All tests passed!"
   end if

contains

   subroutine test_parse_minimal(error)
      !! Test parsing minimal .mqc file with required sections only
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
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
      type(error_t) :: parse_error
      integer :: unit
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
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
      type(error_t) :: parse_error
      integer :: unit
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
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
      type(error_t) :: parse_error
      integer :: unit
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
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
      type(error_t) :: parse_error
      integer :: unit
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
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
      type(error_t) :: parse_error
      integer :: unit
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, allocated(config%log_level), "log_level should be allocated")
      if (allocated(error)) return

      call check(error, config%log_level, "verbose", "log_level should be 'verbose'")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_system_log_level

   subroutine test_parse_with_comments(error)
      !! Test parsing .mqc file with comments after 'end' keywords
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_comments.mqc"

      ! Create test file with comments after 'end' keywords
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end  ! schema'
      write (unit, '(A)') ''
      write (unit, '(A)') '%model'
      write (unit, '(A)') 'method = XTB-GFN2'
      write (unit, '(A)') 'end  ! model'
      write (unit, '(A)') ''
      write (unit, '(A)') '%driver'
      write (unit, '(A)') 'type = Energy'
      write (unit, '(A)') 'end  ! driver'
      write (unit, '(A)') ''
      write (unit, '(A)') '%structure'
      write (unit, '(A)') 'charge = 0'
      write (unit, '(A)') 'multiplicity = 1'
      write (unit, '(A)') 'end  ! structure'
      write (unit, '(A)') ''
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '3'
      write (unit, '(A)') 'Water molecule'
      write (unit, '(A)') 'O 0.0 0.0 0.0'
      write (unit, '(A)') 'H 1.0 0.0 0.0'
      write (unit, '(A)') 'H 0.0 1.0 0.0'
      write (unit, '(A)') 'end  ! geometry'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragments'
      write (unit, '(A)') 'nfrag = 1'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragment'
      write (unit, '(A)') 'charge = 0'
      write (unit, '(A)') 'multiplicity = 1'
      write (unit, '(A)') '%indices'
      write (unit, '(A)') '0 1 2'
      write (unit, '(A)') 'end  ! indices'
      write (unit, '(A)') 'end  ! fragment'
      write (unit, '(A)') ''
      write (unit, '(A)') 'end  ! fragments'
      write (unit, '(A)') ''
      write (unit, '(A)') '%connectivity'
      write (unit, '(A)') 'nbonds = 2'
      write (unit, '(A)') ''
      write (unit, '(A)') '0 1 1 preserved'
      write (unit, '(A)') '0 2 1 preserved'
      write (unit, '(A)') ''
      write (unit, '(A)') 'nbroken = 0'
      write (unit, '(A)') 'end  ! connectivity'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed with comments")
      if (allocated(error)) return

      call check(error, config%schema_name, "mqc-frag", "schema_name should be 'mqc-frag'")
      if (allocated(error)) return

      call check(error, config%method, METHOD_TYPE_GFN2, "method should be GFN2")
      if (allocated(error)) return

      call check(error, config%geometry%natoms, 3, "natoms should be 3")
      if (allocated(error)) return

      call check(error, config%nfrag, 1, "nfrag should be 1")
      if (allocated(error)) return

      call check(error, size(config%fragments(1)%indices), 3, "Fragment should have 3 indices")
      if (allocated(error)) return

      call check(error, config%nbonds, 2, "nbonds should be 2")
      if (allocated(error)) return

      call check(error, config%nbroken, 0, "nbroken should be 0")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_with_comments

   subroutine test_parse_hessian_section(error)
      !! Test parsing %hessian section with custom values
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_hessian.mqc"

      ! Create test file with %hessian section
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%driver'
      write (unit, '(A)') 'type = Hessian'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%hessian'
      write (unit, '(A)') 'displacement = 0.005'
      write (unit, '(A)') 'temperature = 300.0'
      write (unit, '(A)') 'pressure = 1.5'
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%calc_type, CALC_TYPE_HESSIAN, "calc_type should be HESSIAN")
      if (allocated(error)) return

      call check(error, is_equal(config%hessian_displacement, 0.005_dp), &
                 "hessian_displacement should be 0.005")
      if (allocated(error)) return

      call check(error, is_equal(config%hessian_temperature, 300.0_dp), &
                 "hessian_temperature should be 300.0")
      if (allocated(error)) return

      call check(error, is_equal(config%hessian_pressure, 1.5_dp), &
                 "hessian_pressure should be 1.5")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_hessian_section

   subroutine test_parse_hessian_defaults(error)
      !! Test that %hessian section uses defaults when not specified
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_hessian_defaults.mqc"

      ! Create test file without %hessian section
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, is_equal(config%hessian_displacement, DEFAULT_DISPLACEMENT), &
                 "hessian_displacement should be default value")
      if (allocated(error)) return

      call check(error, is_equal(config%hessian_temperature, DEFAULT_TEMPERATURE), &
                 "hessian_temperature should be default value")
      if (allocated(error)) return

      call check(error, is_equal(config%hessian_pressure, DEFAULT_PRESSURE), &
                 "hessian_pressure should be default value")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_hessian_defaults

   subroutine test_parse_aimd_section(error)
      !! Test parsing %aimd section
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_aimd.mqc"

      ! Create test file with %aimd section
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%aimd'
      write (unit, '(A)') 'dt = 0.5'
      write (unit, '(A)') 'nsteps = 1000'
      write (unit, '(A)') 'initial_temperature = 350.0'
      write (unit, '(A)') 'output_frequency = 10'
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, is_equal(config%aimd_dt, 0.5_dp), &
                 "aimd_dt should be 0.5")
      if (allocated(error)) return

      call check(error, config%aimd_nsteps, 1000, "aimd_nsteps should be 1000")
      if (allocated(error)) return

      call check(error, is_equal(config%aimd_initial_temperature, 350.0_dp), &
                 "aimd_initial_temperature should be 350.0")
      if (allocated(error)) return

      call check(error, config%aimd_output_frequency, 10, "aimd_output_frequency should be 10")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_aimd_section

   subroutine test_parse_fragmentation_section(error)
      !! Test parsing %fragmentation section
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_fragmentation.mqc"

      ! Create test file with %fragmentation section
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragmentation'
      write (unit, '(A)') 'method = MBE'
      write (unit, '(A)') 'level = 3'
      write (unit, '(A)') 'allow_overlapping_fragments = true'
      write (unit, '(A)') 'max_intersection_level = 5'
      write (unit, '(A)') 'embedding = none'
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%frag_method, "MBE", "frag_method should be 'MBE'")
      if (allocated(error)) return

      call check(error, config%frag_level, 3, "frag_level should be 3")
      if (allocated(error)) return

      call check(error, config%allow_overlapping_fragments, .true., &
                 "allow_overlapping_fragments should be true")
      if (allocated(error)) return

      call check(error, config%max_intersection_level, 5, "max_intersection_level should be 5")
      if (allocated(error)) return

      call check(error, config%embedding, "none", "embedding should be 'none'")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_fragmentation_section

   subroutine test_parse_fragmentation_cutoffs(error)
      !! Test parsing %fragmentation section with %cutoffs subsection
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_cutoffs.mqc"

      ! Create test file with %fragmentation and %cutoffs
      open (newunit=unit, file=test_file, status='replace', action='write')
      write (unit, '(A)') '%schema'
      write (unit, '(A)') 'name = mqc-frag'
      write (unit, '(A)') 'version = 1.0'
      write (unit, '(A)') 'index_base = 0'
      write (unit, '(A)') 'units = angstrom'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragmentation'
      write (unit, '(A)') 'level = 3'
      write (unit, '(A)') '%cutoffs'
      write (unit, '(A)') '2 = 10.0'
      write (unit, '(A)') '3 = 8.0'
      write (unit, '(A)') 'end'
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%frag_level, 3, "frag_level should be 3")
      if (allocated(error)) return

      call check(error, allocated(config%fragment_cutoffs), "fragment_cutoffs should be allocated")
      if (allocated(error)) return

      call check(error, is_equal(config%fragment_cutoffs(2), 10.0_dp), &
                 "dimer cutoff should be 10.0")
      if (allocated(error)) return

      call check(error, is_equal(config%fragment_cutoffs(3), 8.0_dp), &
                 "trimer cutoff should be 8.0")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_fragmentation_cutoffs

   subroutine test_parse_xtb_solvation(error)
      !! Test parsing %xtb section with solvation settings
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_xtb_solvation.mqc"

      ! Create test file with %xtb solvation settings
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
      write (unit, '(A)') '%xtb'
      write (unit, '(A)') 'solvent = water'
      write (unit, '(A)') 'solvation_model = cpcm'
      write (unit, '(A)') 'use_cds = false'
      write (unit, '(A)') 'use_shift = false'
      write (unit, '(A)') 'dielectric = 80.0'
      write (unit, '(A)') 'cpcm_nang = 230'
      write (unit, '(A)') 'cpcm_rscale = 1.2'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%geometry'
      write (unit, '(A)') '3'
      write (unit, '(A)') 'Water'
      write (unit, '(A)') 'O 0.0 0.0 0.0'
      write (unit, '(A)') 'H 1.0 0.0 0.0'
      write (unit, '(A)') 'H 0.0 1.0 0.0'
      write (unit, '(A)') 'end'
      close (unit)

      ! Parse file
      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%solvent, "water", "solvent should be 'water'")
      if (allocated(error)) return

      call check(error, config%solvation_model, "cpcm", "solvation_model should be 'cpcm'")
      if (allocated(error)) return

      call check(error, config%use_cds, .false., "use_cds should be false")
      if (allocated(error)) return

      call check(error, config%use_shift, .false., "use_shift should be false")
      if (allocated(error)) return

      call check(error, is_equal(config%dielectric, 80.0_dp), &
                 "dielectric should be 80.0")
      if (allocated(error)) return

      call check(error, config%cpcm_nang, 230, "cpcm_nang should be 230")
      if (allocated(error)) return

      call check(error, is_equal(config%cpcm_rscale, 1.2_dp), &
                 "cpcm_rscale should be 1.2")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_parse_xtb_solvation

   subroutine test_parse_fragmentation_global_groups(error)
      !! Test parsing global_groups in %fragmentation section
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_frag_global_groups.mqc"

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
      write (unit, '(A)') '1'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragmentation'
      write (unit, '(A)') 'method = MBE'
      write (unit, '(A)') 'level = 2'
      write (unit, '(A)') 'global_groups = 3'
      write (unit, '(A)') 'end'
      close (unit)

      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%global_groups, 3, "global_groups should be 3")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')
   end subroutine test_parse_fragmentation_global_groups

   subroutine test_parse_fragmentation_nodes_per_group(error)
      !! Test parsing nodes_per_group in %fragmentation section
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_frag_nodes_per_group.mqc"

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
      write (unit, '(A)') '1'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragmentation'
      write (unit, '(A)') 'method = MBE'
      write (unit, '(A)') 'level = 2'
      write (unit, '(A)') 'nodes_per_group = 4'
      write (unit, '(A)') 'end'
      close (unit)

      call read_mqc_file(test_file, config, parse_error)

      call check(error,.not. parse_error%has_error(), "Parser should succeed")
      if (allocated(error)) return

      call check(error, config%nodes_per_group, 4, "nodes_per_group should be 4")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')
   end subroutine test_parse_fragmentation_nodes_per_group

   subroutine test_error_invalid_global_groups(error)
      !! Test invalid global_groups value is rejected
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
      character(len=*), parameter :: test_file = "test_frag_invalid_global_groups.mqc"

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
      write (unit, '(A)') '1'
      write (unit, '(A)') ''
      write (unit, '(A)') 'H 0.0 0.0 0.0'
      write (unit, '(A)') 'end'
      write (unit, '(A)') ''
      write (unit, '(A)') '%fragmentation'
      write (unit, '(A)') 'method = MBE'
      write (unit, '(A)') 'level = 2'
      write (unit, '(A)') 'global_groups = 0'
      write (unit, '(A)') 'end'
      close (unit)

      call read_mqc_file(test_file, config, parse_error)

      call check(error, parse_error%has_error(), "Invalid global_groups should error")
      if (allocated(error)) return

      call check(error, index(parse_error%get_message(), "global_groups") > 0, &
                 "Error message should mention global_groups")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')
   end subroutine test_error_invalid_global_groups

   subroutine test_error_missing_schema(error)
      !! Test error handling when schema section is missing
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error, parse_error%has_error(), "Parser should fail for missing schema")
      if (allocated(error)) return

      call config%destroy()
      open (newunit=unit, file=test_file, status='old', action='read')
      close (unit, status='delete')

   end subroutine test_error_missing_schema

   subroutine test_error_missing_geometry(error)
      !! Test error handling when geometry section is missing
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit
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
      call read_mqc_file(test_file, config, parse_error)

      call check(error, parse_error%has_error(), "Parser should fail for missing geometry")
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
