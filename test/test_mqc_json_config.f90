module test_mqc_json_config
   !! Reads every shipped validation deck and checks what came out.
   !!
   !! This began life as an equivalence test: it read each deck as JSON and as
   !! the `.mqc` that `mqc_prep.py` generated from it, and required the two
   !! configs to agree. That was the argument for trusting the JSON reader
   !! while `.mqc` was still the reference, and it earned its keep -- it caught
   !! a shared local carrying one key's value into another, and a cutoff array
   !! filled with zeros where the reference used a negative sentinel.
   !!
   !! With `.mqc` retired there is no second reader to diff against, so the
   !! expectations are written out: atom counts, fragment counts, bond counts,
   !! charges and multiplicities taken from the decks themselves. Duller than a
   !! cross-check, but they catch the same thing -- a reader that silently
   !! stops seeing part of a real input.
   !!
   !! `test_mqc_json_reader` covers the reader's behaviour on constructed
   !! decks. This one covers it on the real ones, which exercise xyz
   !! resolution, overlapping fragments, connectivity and multi-molecule mode
   !! together, in the combinations users actually write.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use mqc_json_config_reader, only: read_json_config_file
   use mqc_config_types, only: mqc_config_t
   use mqc_error, only: error_t
   use pic_types, only: dp
   implicit none
   private
   public :: collect_mqc_json_config_tests

   character(len=*), parameter :: INPUT_DIR = "../validation/inputs/"
      !! Tests run with their source directory as the working directory

contains

   subroutine collect_mqc_json_config_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("single_molecule_xyz_reference", test_h3o), &
                  new_unittest("fragmented_mbe_with_cutoffs", test_frag_water20), &
                  new_unittest("overlapping_fragments_and_connectivity", test_overlapping_gly3), &
                  new_unittest("multi_molecule_mode", test_multi_frag), &
                  new_unittest("dft_functional", test_dft_b3lyp), &
                  new_unittest("open_shell", test_uhf_oh), &
                  new_unittest("xtb_solvation", test_solvation), &
                  new_unittest("charged_fragments", test_charged_cluster), &
                  new_unittest("fragment_indices_stay_in_range", test_indices_in_range), &
                  new_unittest("missing_file_is_an_error", test_missing_file) &
                  ]
   end subroutine collect_mqc_json_config_tests

   subroutine test_h3o(error)
      !! Geometry reached by xyz path, relative to the deck's own directory
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config

      call load("h3o", config, error)
      if (allocated(error)) return
      call expect_molecule(error, "h3o", config, 0, 4, 0, 0, 1, 1)
      if (allocated(error)) return
      ! The xyz gave real coordinates, not a zeroed array.
      call check(error, maxval(abs(config%geometry%coords)) > 0.0_dp, &
                 "h3o: coordinates should have been read from the xyz file")
   end subroutine test_h3o

   subroutine test_frag_water20(error)
      !! Twenty monomers plus a distance cutoff
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config

      call load("frag_water20_hf", config, error)
      if (allocated(error)) return
      call expect_molecule(error, "frag_water20_hf", config, 0, 60, 20, 0, 0, 1)
      if (allocated(error)) return
      call check(error, allocated(config%fragment_cutoffs), &
                 "frag_water20_hf: cutoffs should be allocated")
      if (allocated(error)) return
      call check(error, abs(config%fragment_cutoffs(2) - 6.0_dp) < 1.0e-10_dp, &
                 "frag_water20_hf: the dimer cutoff should be 6.0")
   end subroutine test_frag_water20

   subroutine test_overlapping_gly3(error)
      !! Overlapping fragments with connectivity, so bonds get broken
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config

      call load("overlapping_gly3", config, error)
      if (allocated(error)) return
      call expect_molecule(error, "overlapping_gly3", config, 0, 24, 3, 3, 0, 1)
      if (allocated(error)) return
      ! These bonds are listed precisely because fragmentation severs them; a
      ! reader that failed to derive is_broken would report none.
      call check(error, config%nbroken > 0, &
                 "overlapping_gly3: fragmentation should break at least one bond")
   end subroutine test_overlapping_gly3

   subroutine test_multi_frag(error)
      !! Two molecules, so the config uses the molecules array
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config

      call load("multi_frag", config, error)
      if (allocated(error)) return
      call check(error, config%nmol, 2, "multi_frag: should be two molecules")
      if (allocated(error)) return
      call check(error, config%molecules(1)%geometry%natoms, 18)
      if (allocated(error)) return
      call check(error, config%molecules(2)%geometry%natoms, 18)
      if (allocated(error)) return
      call check(error, config%molecules(1)%nfrag, 6)
      if (allocated(error)) return
      call check(error, config%molecules(2)%nfrag, 6)
      if (allocated(error)) return
      ! Multi-molecule mode must not also fill the single-molecule fields, or a
      ! consumer keying off nmol would see two molecules and a stray third.
      call check(error, config%geometry%natoms, 0, &
                 "multi_frag: the top-level geometry should stay empty")
   end subroutine test_multi_frag

   subroutine test_dft_b3lyp(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config

      call load("dft_water_b3lyp", config, error)
      if (allocated(error)) return
      call expect_molecule(error, "dft_water_b3lyp", config, 0, 3, 0, 0, 0, 1)
      if (allocated(error)) return
      call check(error, allocated(config%functional), &
                 "dft_water_b3lyp: functional should be set")
      if (allocated(error)) return
      call check(error, allocated(config%basis), "dft_water_b3lyp: basis should be set")
   end subroutine test_dft_b3lyp

   subroutine test_uhf_oh(error)
      !! A doublet: the multiplicity has to survive
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config

      call load("uhf_oh", config, error)
      if (allocated(error)) return
      call expect_molecule(error, "uhf_oh", config, 0, 2, 0, 0, 0, 2)
   end subroutine test_uhf_oh

   subroutine test_solvation(error)
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config

      call load("w1_water_cpcm", config, error)
      if (allocated(error)) return
      call expect_molecule(error, "w1_water_cpcm", config, 0, 4, 0, 0, 1, 1)
      if (allocated(error)) return
      call check(error, allocated(config%solvent), "w1_water_cpcm: solvent should be set")
      if (allocated(error)) return
      call check(error, allocated(config%solvation_model), &
                 "w1_water_cpcm: solvation_model should be set")
   end subroutine test_solvation

   subroutine test_charged_cluster(error)
      !! Fragment charges summing to a non-zero molecular charge
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config
      integer :: ifrag, total

      call load("charged_cluster", config, error)
      if (allocated(error)) return
      call expect_molecule(error, "charged_cluster", config, 0, 28, 8, 0, 4, 1)
      if (allocated(error)) return

      total = 0
      do ifrag = 1, config%nfrag
         total = total + config%fragments(ifrag)%charge
      end do
      call check(error, total, config%charge, &
                 "charged_cluster: fragment charges should sum to the molecular charge")
   end subroutine test_charged_cluster

   subroutine test_indices_in_range(error)
      !! Every fragment index addresses a real atom
      !!
      !! An off-by-one or a dropped index would otherwise surface much later,
      !! as a fragment quietly missing an atom.
      type(error_type), allocatable, intent(out) :: error
      type(mqc_config_t) :: config

      call load("charged_cluster", config, error)
      if (allocated(error)) return
      call expect_indices_in_range(error, "charged_cluster", config)
      if (allocated(error)) return

      call load("overlapping_gly3", config, error)
      if (allocated(error)) return
      call expect_indices_in_range(error, "overlapping_gly3", config)
   end subroutine test_indices_in_range

   subroutine test_missing_file(error)
      !! A path that does not exist must be an error, not an empty config
      type(error_type), allocatable, intent(out) :: error

      type(mqc_config_t) :: config
      type(error_t) :: read_error

      call read_json_config_file(INPUT_DIR//"no_such_deck.json", config, read_error)
      call check(error, read_error%has_error(), &
                 "reading a nonexistent JSON file should set an error")
   end subroutine test_missing_file

   ! ---- helpers -------------------------------------------------------------

   subroutine load(stem, config, error)
      !! Read a shipped deck, failing the test if it does not parse
      character(len=*), intent(in) :: stem
      type(mqc_config_t), intent(out) :: config
      type(error_type), allocatable, intent(out) :: error

      type(error_t) :: read_error

      call read_json_config_file(INPUT_DIR//stem//".json", config, read_error)
      call check(error,.not. read_error%has_error(), &
                 stem//".json failed to parse: "//read_error%get_message())
   end subroutine load

   subroutine expect_molecule(error, stem, config, nmol, natoms, nfrag, nbonds, &
                              charge, multiplicity)
      !! The counts a single-molecule deck should have produced
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: stem
      type(mqc_config_t), intent(in) :: config
      integer, intent(in) :: nmol, natoms, nfrag, nbonds, charge, multiplicity

      call check(error, config%nmol, nmol, stem//": nmol")
      if (allocated(error)) return
      call check(error, config%geometry%natoms, natoms, stem//": natoms")
      if (allocated(error)) return
      call check(error, size(config%geometry%elements), natoms, stem//": element count")
      if (allocated(error)) return
      call check(error, size(config%geometry%coords, 2), natoms, stem//": coordinate count")
      if (allocated(error)) return
      call check(error, config%nfrag, nfrag, stem//": nfrag")
      if (allocated(error)) return
      call check(error, config%nbonds, nbonds, stem//": nbonds")
      if (allocated(error)) return
      call check(error, config%charge, charge, stem//": charge")
      if (allocated(error)) return
      call check(error, config%multiplicity, multiplicity, stem//": multiplicity")
   end subroutine expect_molecule

   subroutine expect_indices_in_range(error, stem, config)
      !! No fragment may name an atom the geometry does not have
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: stem
      type(mqc_config_t), intent(in) :: config

      integer :: ifrag

      do ifrag = 1, config%nfrag
         call check(error, minval(config%fragments(ifrag)%indices) >= 0 .and. &
                    maxval(config%fragments(ifrag)%indices) < config%geometry%natoms, &
                    stem//": fragment index out of range")
         if (allocated(error)) return
      end do
   end subroutine expect_indices_in_range

end module test_mqc_json_config

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_json_config, only: collect_mqc_json_config_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_json_config", collect_mqc_json_config_tests) &
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
