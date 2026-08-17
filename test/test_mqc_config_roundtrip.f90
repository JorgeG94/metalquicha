module test_mqc_config_roundtrip
   !! Guards the config plumbing: input file -> mqc_config_t -> method_config_t
   !! -> concrete method options.
   !!
   !! Every setting is copied by hand at each of those hops, and nothing in the
   !! type system notices when a hop is missed. That has silently dropped
   !! settings three times: `basis` and `aux_basis` were parsed and then
   !! discarded by the adapter, and `functional` was not a recognized keyword
   !! at all -- so an input asking for def2-SVP/PBE quietly ran STO-3G with the
   !! default functional.
   !!
   !! The point of these tests is that a value written in an input file must
   !! arrive, unchanged, at the object that acts on it. Every value below is
   !! deliberately NOT the default, so a dropped copy shows up as the default
   !! rather than coincidentally matching.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_config_types, only: mqc_config_t
   use mqc_json_config_reader, only: read_json_config_file
   use mqc_config_adapter, only: driver_config_t, config_to_driver
   use mqc_method_base, only: qc_method_t
   use mqc_method_factory, only: create_method
   use mqc_method_hf, only: hf_method_t
   use mqc_method_dft, only: dft_method_t
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_config_roundtrip_tests

   character(len=*), parameter :: SCRATCH_FILE = "test_roundtrip_scratch.json"

contains

   subroutine collect_mqc_config_roundtrip_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("counterpoise_roundtrip", test_counterpoise_roundtrip), &
                  new_unittest("hf_settings_reach_the_method", test_hf_roundtrip), &
                  new_unittest("dft_settings_reach_the_method", test_dft_roundtrip) &
                  ]
   end subroutine collect_mqc_config_roundtrip_tests

   subroutine write_input(method, extra_model)
      !! Write a minimal deck with deliberately non-default settings
      character(len=*), intent(in) :: method
      character(len=*), intent(in) :: extra_model  !! e.g. '"functional": "pbe0"'
      integer :: unit

      open (newunit=unit, file=SCRATCH_FILE, status="replace", action="write")
      write (unit, "(A)") "{"
      write (unit, "(A)") '  "schema": {"name": "roundtrip", "version": "1.0"},'
      write (unit, "(A)") '  "model": {'
      write (unit, "(A)") '    "method": "'//method//'",'
      write (unit, "(A)") '    "basis": "cc-pvdz",'
      if (len_trim(extra_model) > 0) then
         write (unit, "(A)") '    "aux_basis": "def2-universal-jkfit",'
         write (unit, "(A)") "    "//extra_model
      else
         write (unit, "(A)") '    "aux_basis": "def2-universal-jkfit"'
      end if
      write (unit, "(A)") "  },"
      write (unit, "(A)") '  "driver": "Energy",'
      write (unit, "(A)") '  "molecules": [{'
      write (unit, "(A)") '    "symbols": ["He"], "geometry": [0.0, 0.0, 0.0],'
      write (unit, "(A)") '    "molecular_charge": 0, "molecular_multiplicity": 1'
      write (unit, "(A)") "  }]"
      write (unit, "(A)") "}"
      close (unit)
   end subroutine write_input

   subroutine write_counterpoise_input(value)
      !! A deck naming a counterpoise scheme, and nothing else unusual
      character(len=*), intent(in) :: value
      integer :: unit

      open (newunit=unit, file=SCRATCH_FILE, status="replace", action="write")
      write (unit, "(A)") "{"
      write (unit, "(A)") '  "schema": {"name": "roundtrip", "version": "1.0"},'
      write (unit, "(A)") '  "model": {"method": "hf", "basis": "cc-pvdz"},'
      write (unit, "(A)") '  "driver": "Energy",'
      write (unit, "(A)") '  "keywords": {"fragmentation": {'
      write (unit, "(A)") '    "method": "MBE", "level": 2,'
      write (unit, "(A)") '    "counterpoise": "'//value//'"'
      write (unit, "(A)") "  }},"
      write (unit, "(A)") '  "molecules": [{'
      write (unit, "(A)") '    "symbols": ["He"], "geometry": [0.0, 0.0, 0.0],'
      write (unit, "(A)") '    "molecular_charge": 0, "molecular_multiplicity": 1'
      write (unit, "(A)") "  }]"
      write (unit, "(A)") "}"
      close (unit)
   end subroutine write_counterpoise_input

   subroutine remove_input()
      integer :: unit
      logical :: exists

      inquire (file=SCRATCH_FILE, exist=exists)
      if (.not. exists) return
      open (newunit=unit, file=SCRATCH_FILE, status="old")
      close (unit, status="delete")
   end subroutine remove_input

   subroutine test_counterpoise_roundtrip(error)
      !! `counterpoise` survives the deck, the config and the adapter
      !!
      !! Three hops, and the middle one is where a keyword usually dies: the
      !! schema validator rejects anything it does not know, so a field added to
      !! the config type but not to the allow-list is refused before the reader
      !! ever sees it.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver
      type(error_t) :: parse_error

      call write_counterpoise_input("vmfc")
      call read_json_config_file(SCRATCH_FILE, config, parse_error)
      call remove_input()

      call check(error,.not. parse_error%has_error(), &
                 "a deck naming counterpoise should parse: "//parse_error%get_full_trace())
      if (allocated(error)) return

      call check(error, allocated(config%counterpoise), &
                 "counterpoise must survive parsing")
      if (allocated(error)) return
      call check(error, trim(config%counterpoise), "vmfc", &
                 "the counterpoise value must survive parsing")
      if (allocated(error)) return

      call config_to_driver(config, driver)
      call check(error, trim(driver%counterpoise), "vmfc", &
                 "counterpoise must survive the config adapter")
      if (allocated(error)) return

      ! And the default is the uncorrected expansion, not an error.
      call write_counterpoise_input("none")
      call read_json_config_file(SCRATCH_FILE, config, parse_error)
      call remove_input()
      call config_to_driver(config, driver)
      call check(error, trim(driver%counterpoise), "none", &
                 "a scheme of none should reach the driver as itself")
   end subroutine test_counterpoise_roundtrip

   subroutine test_hf_roundtrip(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver
      class(qc_method_t), allocatable :: method
      type(error_t) :: parse_error

      call write_input("hf", "")
      call read_json_config_file(SCRATCH_FILE, config, parse_error)
      call remove_input()

      call check(error,.not. parse_error%has_error(), "input should parse")
      if (allocated(error)) return

      ! Hop 1: parser -> mqc_config_t
      call check(error, allocated(config%basis), "basis must survive parsing")
      if (allocated(error)) return
      call check(error, trim(config%basis), "cc-pvdz", "basis value must survive parsing")
      if (allocated(error)) return

      ! Hop 2: mqc_config_t -> method_config_t
      call config_to_driver(config, driver)
      call check(error, trim(driver%method_config%basis_set), "cc-pvdz", &
                 "basis must survive the config adapter")
      if (allocated(error)) return
      call check(error, trim(driver%method_config%scf%aux_basis_set), "def2-universal-jkfit", &
                 "aux_basis must survive the config adapter")
      if (allocated(error)) return

      ! Hop 3: method_config_t -> concrete method options
      allocate (method, source=create_method(driver%method_config))
      select type (m => method)
      type is (hf_method_t)
         call check(error, trim(m%options%basis_set), "cc-pvdz", &
                    "basis must reach the HF method")
         if (allocated(error)) return
         call check(error, trim(m%options%aux_basis_set), "def2-universal-jkfit", &
                    "aux_basis must reach the HF method")
      class default
         call check(error, .false., "method = hf should build an hf_method_t")
      end select
   end subroutine test_hf_roundtrip

   subroutine test_dft_roundtrip(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(driver_config_t) :: driver
      class(qc_method_t), allocatable :: method
      type(error_t) :: parse_error

      call write_input("dft", '"functional": "pbe0"')
      call read_json_config_file(SCRATCH_FILE, config, parse_error)
      call remove_input()

      call check(error,.not. parse_error%has_error(), "input should parse")
      if (allocated(error)) return

      call check(error, allocated(config%functional), "functional must survive parsing")
      if (allocated(error)) return

      call config_to_driver(config, driver)
      call check(error, trim(driver%method_config%dft%functional), "pbe0", &
                 "functional must survive the config adapter")
      if (allocated(error)) return

      allocate (method, source=create_method(driver%method_config))
      select type (m => method)
      type is (dft_method_t)
         call check(error, trim(m%options%functional), "pbe0", &
                    "functional must reach the DFT method")
         if (allocated(error)) return
         call check(error, trim(m%options%basis_set), "cc-pvdz", &
                    "basis must reach the DFT method")
         if (allocated(error)) return
         ! The DFT-specific auxiliary field is empty by default, so this only
         ! passes if the shared SCF value is used as the fallback.
         call check(error, trim(m%options%aux_basis_set), "def2-universal-jkfit", &
                    "aux_basis must reach the DFT method")
      class default
         call check(error, .false., "method = dft should build a dft_method_t")
      end select
   end subroutine test_dft_roundtrip

end module test_mqc_config_roundtrip

program tester_mqc_config_roundtrip
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_config_roundtrip, only: collect_mqc_config_roundtrip_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_config_roundtrip", collect_mqc_config_roundtrip_tests) &
                ]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if
end program tester_mqc_config_roundtrip
