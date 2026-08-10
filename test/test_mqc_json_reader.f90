module test_mqc_json_reader
   !! Behavioural tests for the JSON input reader.
   !!
   !! Ported from the `.mqc` parser's suite when that format was retired: the
   !! assertions are the same requirements, restated against the format that
   !! replaced it. They are what the reader is actually *for*, as distinct from
   !! `test_mqc_json_config`, which checks it against real decks.
   !!
   !! Each test writes a scratch deck and reads it back. The deck is assembled
   !! by `write_deck` from a handful of optional blocks, so a test says what it
   !! varies and nothing else -- the `.mqc` originals repeated forty lines of
   !! file-writing per case, which buried the one line that mattered.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_config_reader, only: read_json_config_file
   use mqc_config_types, only: mqc_config_t
   use mqc_method_types, only: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
   use mqc_calc_types, only: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   use mqc_calculation_defaults, only: DEFAULT_DISPLACEMENT, DEFAULT_TEMPERATURE, &
                                       DEFAULT_PRESSURE
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_json_reader_tests
   public :: remove_deck  !! Tidy the scratch file away

   character(len=*), parameter :: DECK = "test_json_reader_scratch.json"

contains

   subroutine collect_mqc_json_reader_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("minimal_deck", test_minimal), &
                  new_unittest("fragments", test_fragments), &
                  new_unittest("connectivity_marks_broken_bonds", test_connectivity), &
                  new_unittest("no_fragments", test_no_fragments), &
                  new_unittest("xtb_method_spelling", test_method_xtb), &
                  new_unittest("logger_level", test_log_level), &
                  new_unittest("hessian_settings", test_hessian), &
                  new_unittest("hessian_defaults", test_hessian_defaults), &
                  new_unittest("aimd_settings", test_aimd), &
                  new_unittest("fragmentation_settings", test_fragmentation), &
                  new_unittest("fragmentation_cutoffs", test_cutoffs), &
                  new_unittest("cutoffs_must_decrease", test_cutoffs_monotonic), &
                  new_unittest("global_groups", test_global_groups), &
                  new_unittest("nodes_per_group", test_nodes_per_group), &
                  new_unittest("xtb_solvation", test_solvation), &
                  new_unittest("inline_symbols_and_geometry", test_inline_geometry), &
                  new_unittest("error_missing_schema", test_missing_schema), &
                  new_unittest("error_missing_molecules", test_missing_molecules), &
                  new_unittest("error_malformed_json", test_malformed) &
                  ]
   end subroutine collect_mqc_json_reader_tests

   ! ==========================================================================
   !  Deck construction
   ! ==========================================================================

   subroutine write_deck(model, driver, keywords, system, molecule)
      !! Write a scratch deck, omitting any block given as an empty string
      !!
      !! Every argument is a JSON fragment without its trailing comma; the
      !! commas are placed here so a caller never has to think about them.
      character(len=*), intent(in) :: model     !! Body of "model"
      character(len=*), intent(in) :: driver    !! Value of "driver"
      character(len=*), intent(in) :: keywords  !! Body of "keywords", or ""
      character(len=*), intent(in) :: system    !! Body of "system", or ""
      character(len=*), intent(in) :: molecule  !! Body of molecules[0]

      integer :: unit

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") "{"
      write (unit, "(A)") '  "schema": {"name": "mqc-frag", "version": "1.0"},'
      write (unit, "(A)") '  "model": {'//model//'},'
      write (unit, "(A)") '  "driver": "'//driver//'",'
      if (len_trim(keywords) > 0) write (unit, "(A)") '  "keywords": {'//keywords//'},'
      if (len_trim(system) > 0) write (unit, "(A)") '  "system": {'//system//'},'
      write (unit, "(A)") '  "molecules": [{'//molecule//'}]'
      write (unit, "(A)") "}"
      close (unit)
   end subroutine write_deck

   pure function two_atoms() result(text)
      !! A two-atom hydrogen molecule given inline
      character(len=:), allocatable :: text
      text = '"symbols": ["H", "H"], "geometry": [0.0, 0.0, 0.0, 0.7, 0.0, 0.0], '// &
             '"molecular_charge": 0, "molecular_multiplicity": 1'
   end function two_atoms

   subroutine read_deck(config, parse_error)
      !! Read the scratch deck back
      type(mqc_config_t), intent(out) :: config
      type(error_t), intent(out) :: parse_error

      call read_json_config_file(DECK, config, parse_error)
   end subroutine read_deck

   ! ==========================================================================
   !  Tests
   ! ==========================================================================

   subroutine test_minimal(error)
      !! The smallest deck the reader accepts, and every field it implies
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), &
                 "minimal deck should parse: "//parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%schema_name, "mqc-frag")
      if (allocated(error)) return
      call check(error, config%schema_version, "1.0")
      if (allocated(error)) return
      call check(error, config%index_base, 0, "index_base defaults to 0")
      if (allocated(error)) return
      call check(error, config%units, "angstrom", "units default to angstrom")
      if (allocated(error)) return
      call check(error, config%method, METHOD_TYPE_GFN2)
      if (allocated(error)) return
      call check(error, config%calc_type, CALC_TYPE_ENERGY)
      if (allocated(error)) return
      call check(error, config%charge, 0)
      if (allocated(error)) return
      call check(error, config%multiplicity, 1)
      if (allocated(error)) return
      call check(error, config%nmol, 0, "one molecule uses single-molecule mode")
      if (allocated(error)) return
      call check(error, config%geometry%natoms, 2)
      if (allocated(error)) return
      call check(error, trim(config%geometry%elements(1)), "H")
      if (allocated(error)) return
      call check(error, close_enough(config%geometry%coords(1, 2), 0.7_dp), &
                 "second atom sits at x = 0.7")
   end subroutine test_minimal

   subroutine test_fragments(error)
      !! Fragment definitions and their per-fragment charge and multiplicity
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN1"', "Gradient", "", "", &
                      '"symbols": ["H", "H", "O", "H"], '// &
                      '"geometry": [0,0,0, 0.7,0,0, 2,0,0, 2.7,0,0], '// &
                      '"molecular_charge": 1, "molecular_multiplicity": 2, '// &
                      '"fragments": [[0, 1], [2, 3]], '// &
                      '"fragment_charges": [0, 1], '// &
                      '"fragment_multiplicities": [1, 2]')
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%method, METHOD_TYPE_GFN1)
      if (allocated(error)) return
      call check(error, config%calc_type, CALC_TYPE_GRADIENT)
      if (allocated(error)) return
      call check(error, config%charge, 1)
      if (allocated(error)) return
      call check(error, config%multiplicity, 2)
      if (allocated(error)) return
      call check(error, config%nfrag, 2)
      if (allocated(error)) return
      call check(error, size(config%fragments(1)%indices), 2)
      if (allocated(error)) return
      call check(error, config%fragments(1)%indices(1), 0, "indices stay 0-based")
      if (allocated(error)) return
      call check(error, config%fragments(1)%charge, 0)
      if (allocated(error)) return
      call check(error, config%fragments(2)%charge, 1)
      if (allocated(error)) return
      call check(error, config%fragments(2)%multiplicity, 2)
   end subroutine test_fragments

   subroutine test_connectivity(error)
      !! A bond between two fragments is broken; one inside a fragment is not
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", &
                      '"symbols": ["H", "H", "H", "H"], '// &
                      '"geometry": [0,0,0, 0.7,0,0, 2,0,0, 2.7,0,0], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1, '// &
                      '"fragments": [[0, 1], [2, 3]], '// &
                      '"connectivity": [[0, 1, 1], [1, 2, 1]]')
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%nbonds, 2)
      if (allocated(error)) return
      ! 0-1 lies inside fragment 1; 1-2 crosses from fragment 1 to fragment 2.
      call check(error,.not. config%bonds(1)%is_broken, &
                 "a bond inside one fragment is not broken")
      if (allocated(error)) return
      call check(error, config%bonds(2)%is_broken, &
                 "a bond between two fragments is broken")
      if (allocated(error)) return
      call check(error, config%nbroken, 1)
      if (allocated(error)) return
      call check(error, config%bonds(2)%atom_i, 1)
      if (allocated(error)) return
      call check(error, config%bonds(2)%atom_j, 2)
      if (allocated(error)) return
      call check(error, config%bonds(2)%order, 1)
   end subroutine test_connectivity

   subroutine test_no_fragments(error)
      !! A deck without a fragments key is unfragmented, not an error
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%nfrag, 0)
      if (allocated(error)) return
      call check(error, config%nbonds, 0)
   end subroutine test_no_fragments

   subroutine test_method_xtb(error)
      !! "XTB-GFN1" and "gfn1" mean the same method
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN1"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%method, METHOD_TYPE_GFN1)
      if (allocated(error)) return

      call write_deck('"method": "gfn1"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%method, METHOD_TYPE_GFN1, &
                 "the bare spelling should reach the same method")
   end subroutine test_method_xtb

   subroutine test_log_level(error)
      !! system.logger.level reaches log_level, and defaults when absent
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", &
                      '"logger": {"level": "verbose"}', two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%log_level), "log_level should be set")
      if (allocated(error)) return
      call check(error, config%log_level, "verbose")
      if (allocated(error)) return

      ! Absent, the default stands -- and does not pick up a neighbouring key.
      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%log_level, "info", "log_level defaults to info")
      if (allocated(error)) return
      call check(error, config%fragment_breakdown, "csv", &
                 "fragment_breakdown defaults to csv")
   end subroutine test_log_level

   subroutine test_hessian(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Hessian", &
                      '"hessian": {"finite_difference_displacement": 0.005, '// &
                      '"temperature": 300.0, "pressure": 1.5}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%calc_type, CALC_TYPE_HESSIAN)
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_displacement, 0.005_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_temperature, 300.0_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_pressure, 1.5_dp))
   end subroutine test_hessian

   subroutine test_hessian_defaults(error)
      !! An empty hessian block leaves every default in place
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Hessian", '"hessian": {}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_displacement, DEFAULT_DISPLACEMENT))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_temperature, DEFAULT_TEMPERATURE))
      if (allocated(error)) return
      call check(error, close_enough(config%hessian_pressure, DEFAULT_PRESSURE))
   end subroutine test_hessian_defaults

   subroutine test_aimd(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"aimd": {"dt": 0.5, "nsteps": 1000, '// &
                      '"initial_temperature": 350.0, "output_frequency": 10}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%aimd_dt, 0.5_dp))
      if (allocated(error)) return
      call check(error, config%aimd_nsteps, 1000)
      if (allocated(error)) return
      call check(error, close_enough(config%aimd_initial_temperature, 350.0_dp))
      if (allocated(error)) return
      call check(error, config%aimd_output_frequency, 10)
   end subroutine test_aimd

   subroutine test_fragmentation(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "level": 3, '// &
                      '"allow_overlapping_fragments": true, '// &
                      '"max_intersection_level": 5, "embedding": "none", '// &
                      '"cutoff_method": "distance", "distance_metric": "min"}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%frag_method, "MBE")
      if (allocated(error)) return
      call check(error, config%frag_level, 3)
      if (allocated(error)) return
      call check(error, config%allow_overlapping_fragments, .true.)
      if (allocated(error)) return
      call check(error, config%max_intersection_level, 5)
      if (allocated(error)) return
      call check(error, config%embedding, "none")
      if (allocated(error)) return
      call check(error, config%cutoff_method, "distance")
      if (allocated(error)) return
      call check(error, config%distance_metric, "min")
   end subroutine test_fragmentation

   subroutine test_cutoffs(error)
      !! Named and numeric n-mer keys both land at the right level
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "level": 3, '// &
                      '"cutoffs": {"dimer": 10.0, "trimer": 8.0}}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, allocated(config%fragment_cutoffs), "cutoffs should be allocated")
      if (allocated(error)) return
      call check(error, close_enough(config%fragment_cutoffs(2), 10.0_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%fragment_cutoffs(3), 8.0_dp))
      if (allocated(error)) return
      ! A level nobody gave keeps the sentinel, which downstream reads as
      ! "no cutoff here" -- a zero would screen every tetramer out instead.
      call check(error, config%fragment_cutoffs(4) < 0.0_dp, &
                 "an unspecified level keeps the negative sentinel")
      if (allocated(error)) return

      ! The numeric spelling means the same thing.
      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "level": 3, '// &
                      '"cutoffs": {"2": 10.0, "3": 8.0}}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, close_enough(config%fragment_cutoffs(3), 8.0_dp), &
                 "numeric cutoff keys should reach the same level")
   end subroutine test_cutoffs

   subroutine test_cutoffs_monotonic(error)
      !! A trimer cutoff above the dimer one screens in more trimers than
      !! dimers, which cannot be intended
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "level": 3, '// &
                      '"cutoffs": {"dimer": 4.0, "trimer": 9.0}}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error, parse_error%has_error(), &
                 "cutoffs that increase with level should be rejected")
   end subroutine test_cutoffs_monotonic

   subroutine test_global_groups(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "global_groups": 3}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%global_groups, 3)
   end subroutine test_global_groups

   subroutine test_nodes_per_group(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"fragmentation": {"method": "MBE", "nodes_per_group": 4}', &
                      "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%nodes_per_group, 4)
   end subroutine test_nodes_per_group

   subroutine test_solvation(error)
      !! Every xTB solvation knob, including the two only .mqc could reach
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"xtb": {"solvent": "water", "solvation_model": "cpcm", '// &
                      '"dielectric": 80.0, "cpcm_nang": 230, "cpcm_rscale": 1.2, '// &
                      '"use_cds": false, "use_shift": false}', "", two_atoms())
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%solvent, "water")
      if (allocated(error)) return
      call check(error, config%solvation_model, "cpcm")
      if (allocated(error)) return
      call check(error, close_enough(config%dielectric, 80.0_dp))
      if (allocated(error)) return
      call check(error, config%cpcm_nang, 230)
      if (allocated(error)) return
      call check(error, close_enough(config%cpcm_rscale, 1.2_dp))
      if (allocated(error)) return
      call check(error, config%use_cds, .false., "use_cds should be settable from JSON")
      if (allocated(error)) return
      call check(error, config%use_shift, .false., "use_shift should be settable from JSON")
      if (allocated(error)) return

      ! Both default to on, so a deck that says nothing keeps the terms.
      call write_deck('"method": "XTB-GFN2"', "Energy", &
                      '"xtb": {"solvent": "water"}', "", two_atoms())
      call read_deck(config, parse_error)
      call check(error, config%use_cds, .true., "use_cds defaults to true")
      if (allocated(error)) return
      call check(error, config%use_shift, .true., "use_shift defaults to true")
   end subroutine test_solvation

   subroutine test_inline_geometry(error)
      !! symbols + a flat coordinate list, the alternative to an xyz reference
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error

      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", &
                      '"symbols": ["O", "H", "H"], '// &
                      '"geometry": [0.0, 0.0, 0.117, 0.0, 0.757, -0.469, '// &
                      '0.0, -0.757, -0.469], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1')
      call read_deck(config, parse_error)

      call check(error,.not. parse_error%has_error(), parse_error%get_message())
      if (allocated(error)) return
      call check(error, config%geometry%natoms, 3)
      if (allocated(error)) return
      call check(error, trim(config%geometry%elements(1)), "O")
      if (allocated(error)) return
      call check(error, close_enough(config%geometry%coords(3, 1), 0.117_dp))
      if (allocated(error)) return
      call check(error, close_enough(config%geometry%coords(2, 3), -0.757_dp))
      if (allocated(error)) return

      ! A coordinate list that does not match the symbol count is an error,
      ! not a silently truncated molecule.
      call write_deck('"method": "XTB-GFN2"', "Energy", "", "", &
                      '"symbols": ["O", "H", "H"], "geometry": [0.0, 0.0, 0.1], '// &
                      '"molecular_charge": 0, "molecular_multiplicity": 1')
      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), &
                 "a short geometry list should be rejected")
   end subroutine test_inline_geometry

   subroutine test_missing_schema(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") "{"
      write (unit, "(A)") '  "model": {"method": "XTB-GFN2"},'
      write (unit, "(A)") '  "driver": "Energy",'
      write (unit, "(A)") '  "molecules": [{'//two_atoms()//'}]'
      write (unit, "(A)") "}"
      close (unit)

      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "a deck without schema should fail")
   end subroutine test_missing_schema

   subroutine test_missing_molecules(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") "{"
      write (unit, "(A)") '  "schema": {"name": "mqc-frag", "version": "1.0"},'
      write (unit, "(A)") '  "model": {"method": "XTB-GFN2"},'
      write (unit, "(A)") '  "driver": "Energy"'
      write (unit, "(A)") "}"
      close (unit)

      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "a deck without molecules should fail")
   end subroutine test_missing_molecules

   subroutine test_malformed(error)
      !! Broken JSON is a parse error, not a crash or a half-filled config
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      type(error_t) :: parse_error
      integer :: unit

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") '{"schema": {"name": "mqc-frag",'
      close (unit)

      call read_deck(config, parse_error)
      call check(error, parse_error%has_error(), "malformed JSON should fail cleanly")
   end subroutine test_malformed

   pure function close_enough(a, b) result(same)
      real(dp), intent(in) :: a, b
      logical :: same
      same = abs(a - b) <= 1.0e-10_dp*max(1.0_dp, abs(b))
   end function close_enough

   subroutine remove_deck()
      !! Delete the scratch deck
      !!
      !! Called from the driver rather than from each test, so it runs whether
      !! or not a test returned early. It matters more here than for most
      !! scratch files: one case deliberately writes malformed JSON, and a
      !! leftover copy trips the repository's check-json hook.
      integer :: unit
      logical :: exists

      inquire (file=DECK, exist=exists)
      if (.not. exists) return
      open (newunit=unit, file=DECK, status="old")
      close (unit, status="delete")
   end subroutine remove_deck

end module test_mqc_json_reader

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_json_reader, only: collect_mqc_json_reader_tests, remove_deck
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_json_reader", collect_mqc_json_reader_tests) &
                ]

   do is = 1, size(testsuites)
      write (*, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   call remove_deck()

   if (stat > 0) then
      write (error_unit, "(i0, 1x, a)") stat, "test(s) failed!"
      error stop
   end if

end program tester
