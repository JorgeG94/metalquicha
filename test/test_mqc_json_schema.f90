module test_mqc_json_schema
   !! Tests that a malformed deck is rejected, and says why.
   !!
   !! The failure this guards against is silent. Before the validator, a key
   !! the reader did not recognise looked exactly like a key nobody wrote, so
   !! `"max_iter"` for `"maxiter"` ran the default iteration count and reported
   !! nothing. Every rejection below is a case that used to be accepted and
   !! quietly ignored.
   !!
   !! Each test checks the message as well as the rejection. An error that does
   !! not name the offending key leaves the user to bisect their own input,
   !! which for a fragmented deck of several hundred lines is most of the work.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_config_reader, only: read_json_config_file
   use mqc_config_types, only: mqc_config_t
   use mqc_error, only: error_t
   implicit none
   private
   public :: collect_mqc_json_schema_tests
   public :: remove_deck  !! Tidy the scratch file away

   character(len=*), parameter :: DECK = "test_json_schema_scratch.json"

contains

   subroutine collect_mqc_json_schema_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("valid_deck_is_accepted", test_valid), &
                  new_unittest("unknown_root_key", test_unknown_root_key), &
                  new_unittest("misspelled_nested_key", test_misspelled_nested_key), &
                  new_unittest("unknown_molecule_key", test_unknown_molecule_key), &
                  new_unittest("unknown_cutoff_key", test_unknown_cutoff_key), &
                  new_unittest("missing_required_root_key", test_missing_required), &
                  new_unittest("missing_required_nested_key", test_missing_nested_required), &
                  new_unittest("fragment_charges_must_sum", test_charge_sum), &
                  new_unittest("parallel_lists_must_match", test_list_lengths), &
                  new_unittest("geometry_given_two_ways", test_geometry_both), &
                  new_unittest("geometry_given_no_way", test_geometry_neither), &
                  new_unittest("empty_molecules_array", test_empty_molecules), &
                  new_unittest("ecp_is_an_allowed_model_key", test_ecp_allowed), &
                  new_unittest("ecp_refused_for_xtb", test_ecp_xtb), &
                  new_unittest("ecp_refused_for_mcscf", test_ecp_mcscf), &
                  new_unittest("empty_ecp_is_not_a_request", test_ecp_empty) &
                  ]
   end subroutine collect_mqc_json_schema_tests

   ! ---- the cases ----------------------------------------------------------

   subroutine test_valid(error)
      !! Everything the schema allows, together, must still be accepted
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0", "index_base": 0, "units": "angstrom"},', &
                       '"model": {"method": "XTB-GFN2", "basis": "sto-3g", "aux_basis": "x", "functional": "pbe"},', &
                       '"driver": "Energy",', &
                       '"title": "anything",', &
                       '"system": {"logger": {"level": "info"}, "skip_json_output": false,', &
                       '           "fragment_breakdown": "csv"},', &
                       '"keywords": {', &
                       '  "scf": {"maxiter": 40, "tolerance": 1e-8, "unrestricted": false, "guess": "gwh",', &
                       '          "allow_crap_scf": false},', &
                       '  "hessian": {"displacement": 0.005, "temperature": 300.0, "pressure": 1.0},', &
                       '  "aimd": {"dt": 0.5, "nsteps": 10, "initial_temperature": 300.0,', &
                       '           "output_frequency": 1},', &
                       '  "xtb": {"solvent": "water", "solvation_model": "cpcm", "dielectric": 80.0,', &
                       '          "cpcm_nang": 230, "cpcm_rscale": 1.2, "use_cds": true, "use_shift": true},', &
                       '  "fragmentation": {"method": "MBE", "level": 3, "embedding": "none",', &
                       '                    "allow_overlapping_fragments": false,', &
                       '                    "max_intersection_level": 2, "cutoff_method": "distance",', &
                       '                    "distance_metric": "min", "global_groups": 2,', &
                       '                    "cutoffs": {"dimer": 6.0, "3": 4.0}}},', &
                       '"molecules": [{"symbols": ["H", "H", "H", "H"],', &
                       '  "geometry": [0,0,0, 0.7,0,0, 2,0,0, 2.7,0,0],', &
                       '  "molecular_charge": 1, "molecular_multiplicity": 2,', &
                       '  "fragments": [[0, 1], [2, 3]], "fragment_charges": [1, 0],', &
                       '  "fragment_multiplicities": [2, 1],', &
                       '  "connectivity": [[0, 1, 1]]}]'])
      call expect_accepted(error, "a deck using every allowed key")
   end subroutine test_valid

   subroutine test_unknown_root_key(error)
      type(error_type), allocatable, intent(out) :: error

      call write_minimal('"extras": {},')
      call expect_rejected(error, "extras", "an unknown top-level key")
   end subroutine test_unknown_root_key

   subroutine test_misspelled_nested_key(error)
      !! The case that motivated the validator
      type(error_type), allocatable, intent(out) :: error

      call write_minimal('"keywords": {"scf": {"max_iter": 40}},')
      call expect_rejected(error, "max_iter", "a misspelled scf key")
   end subroutine test_misspelled_nested_key

   subroutine test_unknown_molecule_key(error)
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "XTB-GFN2"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"symbols": ["H", "H"], "geometry": [0,0,0, 0.7,0,0],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1, "charge": 0}]'])
      call expect_rejected(error, "charge", "an unknown key on a molecule")
   end subroutine test_unknown_molecule_key

   subroutine test_unknown_cutoff_key(error)
      !! Cutoff keys are n-mer names or level numbers, nothing else
      type(error_type), allocatable, intent(out) :: error

      call write_minimal('"keywords": {"fragmentation": {"method": "MBE", '// &
                         '"cutoffs": {"quadmer": 5.0}}},')
      call expect_rejected(error, "quadmer", "an invented n-mer name")
   end subroutine test_unknown_cutoff_key

   subroutine test_missing_required(error)
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"symbols": ["H"], "geometry": [0,0,0],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1}]'])
      call expect_rejected(error, "model", "a deck with no model block")
   end subroutine test_missing_required

   subroutine test_missing_nested_required(error)
      type(error_type), allocatable, intent(out) :: error

      call write_minimal('"keywords": {"fragmentation": {"level": 2}},')
      call expect_rejected(error, "method", "fragmentation without a method")
   end subroutine test_missing_nested_required

   subroutine test_charge_sum(error)
      !! The check that catches a plausible wrong answer rather than a crash
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "XTB-GFN2"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"symbols": ["H", "H", "H", "H"],', &
                       '  "geometry": [0,0,0, 0.7,0,0, 2,0,0, 2.7,0,0],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1,', &
                       '  "fragments": [[0, 1], [2, 3]], "fragment_charges": [1, 0]}]'])
      call expect_rejected(error, "sum", "fragment charges that do not add up")
   end subroutine test_charge_sum

   subroutine test_list_lengths(error)
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "XTB-GFN2"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"symbols": ["H", "H", "H", "H"],', &
                       '  "geometry": [0,0,0, 0.7,0,0, 2,0,0, 2.7,0,0],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1,', &
                       '  "fragments": [[0, 1], [2, 3]], "fragment_multiplicities": [1]}]'])
      call expect_rejected(error, "fragment_multiplicities", &
                           "a multiplicity list shorter than the fragment list")
   end subroutine test_list_lengths

   subroutine test_geometry_both(error)
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "XTB-GFN2"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"xyz": "somewhere.xyz", "symbols": ["H"], "geometry": [0,0,0],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1}]'])
      call expect_rejected(error, "not both", "a molecule giving xyz and symbols")
   end subroutine test_geometry_both

   subroutine test_geometry_neither(error)
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "XTB-GFN2"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"molecular_charge": 0, "molecular_multiplicity": 1}]'])
      call expect_rejected(error, "either", "a molecule with no geometry at all")
   end subroutine test_geometry_neither

   subroutine test_empty_molecules(error)
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "XTB-GFN2"},', &
                       '"driver": "Energy",', &
                       '"molecules": []'])
      call expect_rejected(error, "molecules", "an empty molecules array")
   end subroutine test_empty_molecules

   ! ---- scaffolding ---------------------------------------------------------

   subroutine write_deck(body)
      !! Write a deck from its body lines, which supply their own commas
      character(len=*), intent(in) :: body(:)
      integer :: unit, i

      open (newunit=unit, file=DECK, status="replace", action="write")
      write (unit, "(A)") "{"
      do i = 1, size(body)
         write (unit, "(A)") "  "//trim(body(i))
      end do
      write (unit, "(A)") "}"
      close (unit)
   end subroutine write_deck

   subroutine write_minimal(extra)
      !! A deck that would be valid but for `extra`, which is inserted verbatim
      character(len=*), intent(in) :: extra

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "XTB-GFN2"},', &
                       '"driver": "Energy",', &
                       extra, &
                       '"molecules": [{"symbols": ["H", "H"], "geometry": [0,0,0, 0.7,0,0],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1}]'])
   end subroutine write_minimal

   subroutine expect_accepted(error, what)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: what

      type(mqc_config_t) :: config
      type(error_t) :: read_error

      call read_json_config_file(DECK, config, read_error)
      call check(error,.not. read_error%has_error(), &
                 what//" should be accepted, but: "//read_error%get_message())
   end subroutine expect_accepted

   ! ---- effective core potentials ------------------------------------------
   !
   ! `model.ecp` reaches the Hartree-Fock and DFT builders and nothing else, so
   ! a deck naming one for a method that does not carry it would run
   ! all-electron and report a number with nothing to say half the request was
   ! dropped. Wrong by hundreds of Hartree, converged, and silent -- which is
   ! the failure the whole ECP path exists to prevent, so it must not be the one
   ! the deck layer has. Each case below is a refusal rather than an omission.

   subroutine test_ecp_allowed(error)
      !! The key itself validates: the allow-list rejects what it does not know,
      !! so `ecp` being absent from it would refuse every ECP deck ever written.
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "hf", "basis": "def2-svp", "ecp": "def2-ecp"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"symbols": ["I", "H"], "geometry": [0,0,0, 0,0,1.609],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1}]'])
      call expect_accepted(error, "an ECP named on a Hartree-Fock deck")
   end subroutine test_ecp_allowed

   subroutine test_ecp_xtb(error)
      !! A semiempirical method has its own parameterised core and no
      !! all-electron basis for a potential to replace, so naming one means the
      !! deck is confused rather than merely redundant.
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "GFN2", "ecp": "def2-ecp"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"symbols": ["H", "H"], "geometry": [0,0,0, 0.7,0,0],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1}]'])
      call expect_rejected(error, "semiempirical", "an ECP on an xTB deck")
   end subroutine test_ecp_xtb

   subroutine test_ecp_mcscf(error)
      !! The MCSCF path does not carry a potential yet. Refused rather than
      !! ignored, which is the difference between an error and a wrong answer.
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "casscf", "basis": "def2-svp", "ecp": "def2-ecp"},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"symbols": ["I", "H"], "geometry": [0,0,0, 0,0,1.609],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1}]'])
      call expect_rejected(error, "MCSCF", "an ECP on a CASSCF deck")
   end subroutine test_ecp_mcscf

   subroutine test_ecp_empty(error)
      !! An empty string is not a request for a potential, so it must not
      !! refuse a method that has no ECP support. The check reads the value and
      !! not merely the key's presence.
      type(error_type), allocatable, intent(out) :: error

      call write_deck([character(len=200) :: &
                       '"schema": {"name": "t", "version": "1.0"},', &
                       '"model": {"method": "GFN2", "ecp": ""},', &
                       '"driver": "Energy",', &
                       '"molecules": [{"symbols": ["H", "H"], "geometry": [0,0,0, 0.7,0,0],', &
                       '  "molecular_charge": 0, "molecular_multiplicity": 1}]'])
      call expect_accepted(error, "an empty ecp on an xTB deck")
   end subroutine test_ecp_empty

   subroutine expect_rejected(error, needle, what)
      !! Rejected, and the message names the thing that was wrong
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: needle  !! Text the message must contain
      character(len=*), intent(in) :: what

      type(mqc_config_t) :: config
      type(error_t) :: read_error
      character(len=:), allocatable :: message

      call read_json_config_file(DECK, config, read_error)
      call check(error, read_error%has_error(), what//" should be rejected")
      if (allocated(error)) return

      message = read_error%get_message()
      call check(error, index(message, needle) > 0, &
                 what//": the message should mention '"//needle//"' but said: "//message)
   end subroutine expect_rejected

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

end module test_mqc_json_schema

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_json_schema, only: collect_mqc_json_schema_tests, remove_deck
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_json_schema", collect_mqc_json_schema_tests) &
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
