!! The output document, written and read back
module test_mqc_json_writer
   !! `write_json_output` is the only place this program writes JSON, and what
   !! it writes is the contract every consumer reads: the validation harness
   !! scrapes `total_energy` from it, the Python interface reads the gradient
   !! norm, the SAPT and bonding sections out of it, and a restart compares its
   !! fingerprint. Nothing checked that contract -- the writer was reached only
   !! by running whole calculations, which assert on energies rather than on
   !! the shape of the file, so a renamed key or a section written under the
   !! wrong mode broke a consumer and no test.
   !!
   !! Each case here fills a `json_output_data_t` by hand, writes it, and reads
   !! the file back with json-fortran. That is the round trip a consumer makes,
   !! and it is deliberately done through `write_json_output` rather than the
   !! mode-specific routines: which writer a mode dispatches to is part of what
   !! is being tested.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int64
   use mqc_json_output_types, only: json_output_data_t, OUTPUT_MODE_UNFRAGMENTED, &
                                    OUTPUT_MODE_MBE, OUTPUT_MODE_GMBE_PIE
   use mqc_json_writer, only: write_json_output
   use mqc_io_helpers, only: set_output_json_filename, get_output_json_filename
   use json_module, only: json_file
   use mqc_program_limits, only: N_EFMO_TERMS
   use mqc_result_types, only: STATE_SPIN_SINGLET, STATE_SPIN_TRIPLET, &
                               STATE_SPIN_UNRESTRICTED, STATE_SPIN_UNKNOWN
   implicit none
   private

   public :: collect_mqc_json_writer_tests

contains

   subroutine collect_mqc_json_writer_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("unfragmented_document_has_its_keys", test_unfragmented), &
                  new_unittest("gradient_and_hessian_norms_are_written", test_derivatives), &
                  new_unittest("dipole_is_written_with_its_magnitude", test_dipole), &
                  new_unittest("mbe_document_carries_its_levels", test_mbe), &
                  new_unittest("pie_document_counts_nonzero_terms", test_pie), &
                  new_unittest("pie_atom_set_with_no_sentinel_stays_in_bounds", test_pie_full_set), &
                  new_unittest("a_fingerprint_is_written_when_there_is_one", test_fingerprint), &
                  new_unittest("excited_states_round_trip", test_excited_states), &
                  new_unittest("efmo_pair_map_round_trips_strongest_first", &
                               test_efmo_pairs), &
                  new_unittest("unrestricted_roots_carry_their_own_spin_word", &
                               test_unrestricted_spin) &
                  ]
   end subroutine collect_mqc_json_writer_tests

   subroutine written_document(data, json, path)
      !! Write `data` and open what came out, under a name of this test's own.
      type(json_output_data_t), intent(inout) :: data
      type(json_file), intent(out) :: json
      character(len=*), intent(in) :: path

      call set_output_json_filename(path)
      call write_json_output(data)
      call json%initialize()
      call json%load_file(trim(get_output_json_filename()))
   end subroutine written_document

   subroutine test_unfragmented(error)
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      real(dp) :: energy, gap
      logical :: found

      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -76.026760737428_dp
      data%has_energy = .true.
      data%has_orbitals = .true.
      data%homo = -0.4919_dp
      data%lumo = 0.1857_dp

      call written_document(data, json, "jw_unfragmented.json")

      call json%get("jw_unfragmented.total_energy", energy, found)
      call check(error, found, "total_energy is missing from the document")
      if (allocated(error)) return
      call check(error, abs(energy - data%total_energy) < 1.0e-12_dp, &
                 "total_energy came back changed")
      if (allocated(error)) return

      ! The gap is written in eV while the orbitals are in Hartree, and the
      ! conversion is the writer's own -- a consumer reading it as Hartree gets
      ! a number 27 times too large and no complaint.
      call json%get("jw_unfragmented.homo_lumo_gap_ev", gap, found)
      call check(error, found, "homo_lumo_gap_ev is missing")
      if (allocated(error)) return
      call check(error, abs(gap - (data%lumo - data%homo)*27.211386245988_dp) < 1.0e-6_dp, &
                 "the gap is not the orbital difference in eV")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()
   end subroutine test_unfragmented

   subroutine test_derivatives(error)
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      real(dp) :: norm
      logical :: found

      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -1.0_dp
      data%has_energy = .true.

      ! A gradient whose norm is exactly 5: 3-4-0 on one atom, nothing on the
      ! other, so a wrong reduction shows up as a round number that is wrong.
      allocate (data%gradient(3, 2))
      data%gradient = 0.0_dp
      data%gradient(1, 1) = 3.0_dp
      data%gradient(2, 1) = 4.0_dp
      data%has_gradient = .true.

      allocate (data%hessian(2, 2))
      data%hessian = 0.0_dp
      data%hessian(1, 1) = 3.0_dp
      data%hessian(2, 2) = 4.0_dp
      data%has_hessian = .true.

      call written_document(data, json, "jw_derivatives.json")

      call json%get("jw_derivatives.gradient_norm", norm, found)
      call check(error, found, "gradient_norm is missing")
      if (allocated(error)) return
      call check(error, abs(norm - 5.0_dp) < 1.0e-12_dp, "the gradient norm is wrong")
      if (allocated(error)) return

      call json%get("jw_derivatives.hessian_frobenius_norm", norm, found)
      call check(error, found, "hessian_frobenius_norm is missing")
      if (allocated(error)) return
      call check(error, abs(norm - 5.0_dp) < 1.0e-12_dp, "the Hessian norm is wrong")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()
   end subroutine test_derivatives

   subroutine test_dipole(error)
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      real(dp) :: x, magnitude
      logical :: found

      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -1.0_dp
      data%has_energy = .true.
      allocate (data%dipole(3))
      data%dipole = [0.0_dp, 3.0_dp, 4.0_dp]
      data%has_dipole = .true.

      call written_document(data, json, "jw_dipole.json")

      call json%get("jw_dipole.dipole.x", x, found)
      call check(error, found, "the dipole object is missing")
      if (allocated(error)) return
      call check(error, abs(x) < 1.0e-12_dp, "the dipole components are permuted")
      if (allocated(error)) return

      ! Stored in atomic units and reported in Debye. The factor is the writer's
      ! and a consumer cannot tell from the number which one it got.
      call json%get("jw_dipole.dipole.magnitude_debye", magnitude, found)
      call check(error, found, "the dipole magnitude is missing")
      if (allocated(error)) return
      call check(error, abs(magnitude - 5.0_dp*2.541746_dp) < 1.0e-5_dp, &
                 "the dipole magnitude is not the norm converted to Debye")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()
   end subroutine test_dipole

   subroutine test_mbe(error)
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      real(dp) :: total, monomer_sum
      integer :: count
      logical :: found

      ! Two monomers and the dimer they form: the smallest expansion that has
      ! more than one level to lay out.
      data%output_mode = OUTPUT_MODE_MBE
      data%total_energy = -152.0_dp
      data%has_energy = .true.
      data%fragment_count = 3_int64
      data%max_level = 2
      data%fragment_breakdown = "none"   ! no CSV beside the document here
      allocate (data%polymers(3, 2))
      data%polymers = 0
      data%polymers(1, 1) = 1
      data%polymers(2, 1) = 2
      data%polymers(3, 1) = 1
      data%polymers(3, 2) = 2
      allocate (data%fragment_energies(3))
      data%fragment_energies = [-76.0_dp, -76.0_dp, -152.001_dp]
      allocate (data%delta_energies(3))
      data%delta_energies = [0.0_dp, 0.0_dp, -0.001_dp]
      allocate (data%fragment_distances(3))
      data%fragment_distances = [0.0_dp, 0.0_dp, 2.9_dp]
      allocate (data%sum_by_level(2))
      data%sum_by_level = [-152.0_dp, -0.001_dp]

      call written_document(data, json, "jw_mbe.json")

      call json%get("jw_mbe.total_energy", total, found)
      call check(error, found, "total_energy is missing from an MBE document")
      if (allocated(error)) return

      ! The levels array is what a consumer walks, and the monomer level must
      ! carry the count it was given rather than the fragment total.
      call json%get("jw_mbe.levels(1).count", count, found)
      call check(error, found, "the levels array is missing")
      if (allocated(error)) return
      call check(error, count == 2, "the monomer level does not hold two fragments")
      if (allocated(error)) return

      call json%get("jw_mbe.levels(1).total_energy", monomer_sum, found)
      call check(error, found, "a level carries no total")
      if (allocated(error)) return
      call check(error, abs(monomer_sum + 152.0_dp) < 1.0e-9_dp, &
                 "the monomer level total is not the sum it was given")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()
   end subroutine test_mbe

   subroutine test_pie(error)
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      integer :: count
      logical :: found

      ! Three unique atom sets, one of which cancelled to a zero coefficient.
      ! The count is of the terms that survive, which is the number a GMBE
      ! reader compares against its own enumeration.
      data%output_mode = OUTPUT_MODE_GMBE_PIE
      data%total_energy = -228.0_dp
      data%has_energy = .true.
      data%n_pie_terms = 3_int64
      allocate (data%pie_atom_sets(3, 3))
      data%pie_atom_sets = 0
      data%pie_atom_sets(1, 1) = 1
      data%pie_atom_sets(1, 2) = 2
      data%pie_atom_sets(1, 3) = 3
      allocate (data%pie_coefficients(3))
      data%pie_coefficients = [1, 0, -1]
      allocate (data%pie_energies(3))
      data%pie_energies = [-76.0_dp, -76.0_dp, -76.0_dp]

      call written_document(data, json, "jw_pie.json")

      call json%get("jw_pie.pie_terms.count", count, found)
      call check(error, found, "the pie_terms object is missing")
      if (allocated(error)) return
      call check(error, count == 2, &
                 "a term whose coefficient cancelled was counted anyway")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()
   end subroutine test_pie

   subroutine test_pie_full_set(error)
      !! An atom set that fills its column has no negative sentinel to stop on
      !!
      !! The walk that measures each term's atom list is bounded by `max_atoms`
      !! *and* by a negative sentinel. Written as one `.and.` condition, the
      !! bound does not protect the subscript: Fortran may evaluate both
      !! operands, so a column with no sentinel reads `pie_atom_sets(max_atoms
      !! + 1, i)` on the last pass. `-fcheck=bounds` traps it; a release build
      !! reads out of the next column and says nothing.
      !!
      !! Every atom index here is non-negative, so the sentinel never fires and
      !! only the bound can end the walk.
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      integer :: count
      logical :: found

      data%output_mode = OUTPUT_MODE_GMBE_PIE
      data%total_energy = -152.0_dp
      data%has_energy = .true.
      data%n_pie_terms = 1_int64
      allocate (data%pie_atom_sets(2, 1))
      data%pie_atom_sets(:, 1) = [0, 1]
      allocate (data%pie_coefficients(1))
      data%pie_coefficients = [1]
      allocate (data%pie_energies(1))
      data%pie_energies = [-76.0_dp]

      call written_document(data, json, "jw_pie_full.json")

      call json%get("jw_pie_full.pie_terms.count", count, found)
      call check(error, found, "the pie_terms object is missing")
      if (allocated(error)) return
      call check(error, count == 1, "the only term was not written")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()
   end subroutine test_pie_full_set

   subroutine test_fingerprint(error)
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      character(len=:), allocatable :: stamp
      logical :: found

      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -1.0_dp
      data%has_energy = .true.
      data%fingerprint = "0123456789abcdef"

      call written_document(data, json, "jw_fingerprint.json")

      ! What a restart compares before it reuses anything, so an unwritten
      ! fingerprint is a checkpoint that cannot be validated.
      call json%get("jw_fingerprint.fingerprint", stamp, found)
      call check(error, found, "the fingerprint is missing")
      if (allocated(error)) return
      call check(error, stamp == "0123456789abcdef", "the fingerprint came back changed")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()

      ! And absent when there is none: an empty stamp must not be written as a
      ! key holding an empty string, which reads as "checked and matched".
      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -1.0_dp
      data%has_energy = .true.
      data%fingerprint = ""
      call written_document(data, json, "jw_no_fingerprint.json")
      call json%get("jw_no_fingerprint.fingerprint", stamp, found)
      call check(error,.not. found, "an empty fingerprint was written as a key")
      call json%destroy()
      call data%destroy()
   end subroutine test_fingerprint

   subroutine test_excited_states(error)
      !! The spectrum, written and read back state by state
      !!
      !! The per-state object is the contract: a consumer looking for the
      !! brightest root reads one object, not four parallel arrays it has to
      !! index consistently. Both the eV conversion and the spin word are
      !! produced by the writer and exist nowhere in the data, so both are
      !! checked here rather than assumed.
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      real(dp) :: value
      character(len=:), allocatable :: text
      logical :: found

      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -76.026767997_dp
      data%has_energy = .true.
      data%excitation_energies = [0.3386923781_dp, 0.3047529969_dp]
      data%excited_total_energies = [-75.6880756189_dp, -75.7220150001_dp]
      data%oscillator_strengths = [0.02847955_dp, 0.0_dp]
      data%oscillator_strengths_velocity = [0.12950160_dp, 0.0_dp]
      data%transition_dipoles = reshape([0.3551480624_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp], [3, 2])
      data%transition_velocities = reshape([0.2159163096_dp, 0.0_dp, 0.0_dp, &
                                            0.0_dp, 0.0_dp, 0.0_dp], [3, 2])
      data%transition_dipole_origin = [0.0_dp, 0.0_dp, 0.1257326_dp]
      data%nto_leading_weight = [0.9997741843_dp, 0.9812_dp]
      data%state_spin = [STATE_SPIN_SINGLET, STATE_SPIN_TRIPLET]
      data%excited_method = "tda"
      data%excited_spin = "both"
      data%has_excited_states = .true.

      call written_document(data, json, "jw_excited.json")

      call json%get("jw_excited.excited_states.n_states", value, found)
      call check(error, found, "excited_states.n_states is missing")
      if (allocated(error)) return
      call check(error, nint(value) == 2, "n_states is not the number of roots written")
      if (allocated(error)) return

      call json%get("jw_excited.excited_states.method", text, found)
      call check(error, found, "excited_states.method is missing")
      if (allocated(error)) return
      call check(error, text == "tda", "the response problem came back changed")
      if (allocated(error)) return

      call json%get("jw_excited.excited_states.spin", text, found)
      call check(error, found, "excited_states.spin is missing")
      if (allocated(error)) return
      call check(error, text == "both", "the requested spin came back changed")
      if (allocated(error)) return

      call json%get("jw_excited.excited_states.states(1).state", value, found)
      call check(error, found, "the first state object is missing")
      if (allocated(error)) return
      call check(error, nint(value) == 1, "states are numbered from one")
      if (allocated(error)) return

      call json%get("jw_excited.excited_states.states(1).excitation_energy_hartree", &
                    value, found)
      call check(error, found, "excitation_energy_hartree is missing")
      if (allocated(error)) return
      call check(error, abs(value - 0.3386923781_dp) < 1.0e-12_dp, &
                 "the excitation energy came back changed")
      if (allocated(error)) return

      ! Written in eV as well as Hartree, and the conversion is the writer's
      ! own -- a consumer reading the eV column as Hartree is off by 27.
      call json%get("jw_excited.excited_states.states(1).excitation_energy_ev", value, found)
      call check(error, found, "excitation_energy_ev is missing")
      if (allocated(error)) return
      call check(error, abs(value - 0.3386923781_dp*27.211386245988_dp) < 1.0e-6_dp, &
                 "the eV column is not the Hartree one converted")
      if (allocated(error)) return

      ! The state's own total energy, which is what a spectrum is plotted
      ! against the ground state with. It is carried rather than derived:
      ! the reference it sits on is not always the one printed beside it.
      call json%get("jw_excited.excited_states.states(1).total_energy_hartree", &
                    value, found)
      call check(error, found, "total_energy_hartree is missing")
      if (allocated(error)) return
      call check(error, abs(value + 75.6880756189_dp) < 1.0e-12_dp, &
                 "the excited-state total energy came back changed")
      if (allocated(error)) return

      call json%get("jw_excited.excited_states.states(1).spin", text, found)
      call check(error, found, "the state spin label is missing")
      if (allocated(error)) return
      call check(error, text == "singlet", "the first state should be labelled singlet")
      if (allocated(error)) return

      call json%get("jw_excited.excited_states.states(1).oscillator_strength", value, found)
      call check(error, found, "oscillator_strength is missing")
      if (allocated(error)) return
      call check(error, abs(value - 0.02847955_dp) < 1.0e-12_dp, &
                 "the oscillator strength came back changed")
      if (allocated(error)) return

      ! Both gauges, because they are different numbers and a consumer that
      ! read one for the other would be wrong by a factor of four here.
      call json%get("jw_excited.excited_states.states(1).oscillator_strength_velocity", &
                    value, found)
      call check(error, found, "oscillator_strength_velocity is missing")
      if (allocated(error)) return
      call check(error, abs(value - 0.12950160_dp) < 1.0e-12_dp, &
                 "the velocity-gauge oscillator strength came back changed")
      if (allocated(error)) return

      call json%get("jw_excited.excited_states.states(1).transition_dipole(1)", value, found)
      call check(error, found, "the transition dipole is missing")
      if (allocated(error)) return
      call check(error, abs(value - 0.3551480624_dp) < 1.0e-12_dp, &
                 "the transition dipole x component came back changed")
      if (allocated(error)) return

      ! The velocity-gauge moment, beside its length-gauge partner. The two
      ! gauges are the diagnostic the module reports both for, and a document
      ! carrying only one of them cannot show the gap.
      call json%get("jw_excited.excited_states.states(1).transition_velocity(1)", &
                    value, found)
      call check(error, found, "the velocity-gauge transition moment is missing")
      if (allocated(error)) return
      call check(error, abs(value - 0.2159163096_dp) < 1.0e-12_dp, &
                 "the velocity-gauge moment x component came back changed")
      if (allocated(error)) return

      ! One origin for the whole spectrum, on the section rather than per
      ! state: it says which convention the dipoles were measured in.
      call json%get("jw_excited.excited_states.dipole_origin_bohr(3)", value, found)
      call check(error, found, "the transition dipole origin is missing")
      if (allocated(error)) return
      call check(error, abs(value - 0.1257326_dp) < 1.0e-12_dp, &
                 "the transition dipole origin came back changed")
      if (allocated(error)) return

      call json%get("jw_excited.excited_states.states(1).nto_leading_weight", &
                    value, found)
      call check(error, found, "nto_leading_weight is missing")
      if (allocated(error)) return
      call check(error, abs(value - 0.9997741843_dp) < 1.0e-12_dp, &
                 "the leading natural transition orbital weight came back changed")
      if (allocated(error)) return

      ! The second root is a triplet, whose zero oscillator strength is a real
      ! value rather than a missing one.
      call json%get("jw_excited.excited_states.states(2).spin", text, found)
      call check(error, found, "the second state spin label is missing")
      if (allocated(error)) return
      call check(error, text == "triplet", "the second state should be labelled triplet")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()

      ! And nothing at all when no states were computed: the section is the
      ! signal that a spectrum exists.
      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -1.0_dp
      data%has_energy = .true.
      call written_document(data, json, "jw_no_excited.json")
      call json%get("jw_no_excited.excited_states.n_states", value, found)
      call check(error,.not. found, "an excited_states section appeared with no states")
      call json%destroy()
      call data%destroy()
   end subroutine test_excited_states

   subroutine test_efmo_pairs(error)
      !! The per-pair interaction map, written and read back
      !!
      !! Three things exist only in the writer and so are checked rather than
      !! assumed: the ordering, which is by descending magnitude and not the
      !! order the pairs arrive in; the treatment word, which is produced from
      !! a logical; and the suppression of the four named terms on a quantum
      !! pair, which has no such decomposition to report.
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      real(dp) :: value
      character(len=:), allocatable :: text
      logical :: found

      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -228.1_dp
      data%has_energy = .true.
      allocate (data%efmo_terms(N_EFMO_TERMS), source=0.0_dp)
      data%has_efmo = .true.

      ! Deliberately not in magnitude order, and with the largest last, so a
      ! writer that simply echoed the array would fail the first check.
      data%efmo_pair_fragments = reshape([1, 2, 1, 3, 2, 3], [2, 3])
      data%efmo_pair_distance = [0.91_dp, 2.60_dp, 1.75_dp]
      data%efmo_pair_qm = [.true., .false., .false.]
      data%efmo_pair_energy = [-0.0012_dp, -0.0004_dp, -0.0250_dp]
      data%efmo_pair_terms = reshape([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                      -0.0003_dp, -0.0001_dp, 0.0002_dp, -0.0002_dp, &
                                      -0.0200_dp, -0.0030_dp, 0.0010_dp, -0.0030_dp], &
                                     [4, 3])

      call written_document(data, json, "jw_efmo_pairs.json")

      ! Strongest first: the 2-3 pair, written third, has to come back first.
      call json%get("jw_efmo_pairs.efmo.pairs(1).interaction_energy", value, found)
      call check(error, found, "the first pair object is missing")
      if (allocated(error)) return
      call check(error, value, -0.0250_dp, thr=1.0e-12_dp, &
                 message="the pairs did not come back strongest first")
      if (allocated(error)) return

      call json%get("jw_efmo_pairs.efmo.pairs(1).fragments(1)", value, found)
      call check(error, found, "the fragment list is missing")
      if (allocated(error)) return
      call check(error, nint(value) == 2, "the fragment numbers do not follow the sort")
      if (allocated(error)) return

      call json%get("jw_efmo_pairs.efmo.pairs(1).treatment", text, found)
      call check(error, found, "the treatment word is missing")
      if (allocated(error)) return
      call check(error, text == "classical", "a far pair is not called classical")
      if (allocated(error)) return

      call json%get("jw_efmo_pairs.efmo.pairs(1).electrostatics", value, found)
      call check(error, found, "a far pair did not carry its named terms")
      if (allocated(error)) return
      call check(error, value, -0.0200_dp, thr=1.0e-12_dp, &
                 message="the electrostatics term came back changed")
      if (allocated(error)) return

      ! By magnitude the order is 2-3, then the quantum 1-2, then 1-3, so the
      ! quantum pair lands in the middle rather than where it was written.
      call json%get("jw_efmo_pairs.efmo.pairs(2).treatment", text, found)
      call check(error, found, "the middle pair object is missing")
      if (allocated(error)) return
      call check(error, text == "quantum", "the dimer-SCF pair is not called quantum")
      if (allocated(error)) return

      call json%get("jw_efmo_pairs.efmo.pairs(2).electrostatics", value, found)
      call check(error,.not. found, &
                 "a quantum pair reported an electrostatics term it does not have")
      if (allocated(error)) return

      ! And its own distance travelled with it through the sort.
      call json%get("jw_efmo_pairs.efmo.pairs(2).distance", value, found)
      call check(error, found, "the distance is missing")
      if (allocated(error)) return
      call check(error, value, 0.91_dp, thr=1.0e-12_dp, &
                 message="the distance did not travel with its own pair")
      if (allocated(error)) return

      call json%get("jw_efmo_pairs.efmo.pairs(3).interaction_energy", value, found)
      call check(error, found, "the last pair object is missing")
      if (allocated(error)) return
      call check(error, value, -0.0004_dp, thr=1.0e-12_dp, &
                 message="the weakest pair is not last")
   end subroutine test_efmo_pairs

   subroutine test_unrestricted_spin(error)
      !! `STATE_SPIN_UNRESTRICTED` comes back as the word, and so does the gap
      !!
      !! The third spin label, and the one a deck can never ask for: it says
      !! the reference was unrestricted and its roots are not spin
      !! eigenstates. Like the other two it exists only in the writer -- the
      !! container carries an integer -- so a round trip is the only place it
      !! is checked. `STATE_SPIN_UNKNOWN` is here for the same reason: the
      !! `select case` has a default arm whose job is to invent nothing, and
      !! a label silently replaced by a guess would be believed.
      type(error_type), allocatable, intent(out) :: error

      type(json_output_data_t) :: data
      type(json_file) :: json
      character(len=:), allocatable :: text
      logical :: found

      data%output_mode = OUTPUT_MODE_UNFRAGMENTED
      data%total_energy = -75.939126000634_dp
      data%has_energy = .true.
      data%excitation_energies = [0.088470327358_dp, 0.234583947377_dp]
      data%state_spin = [STATE_SPIN_UNRESTRICTED, STATE_SPIN_UNKNOWN]
      data%excited_method = "tda"
      data%excited_spin = "singlet"
      data%has_excited_states = .true.

      call written_document(data, json, "jw_unrestricted.json")

      call json%get("jw_unrestricted.excited_states.states(1).spin", text, found)
      call check(error, found, "the unrestricted state spin label is missing")
      if (allocated(error)) return
      call check(error, text == "unrestricted", "an unrestricted root should be "// &
                 "labelled unrestricted, not with a multiplicity it does not have")
      if (allocated(error)) return

      call json%get("jw_unrestricted.excited_states.states(2).spin", text, found)
      call check(error, found, "the unassigned state spin label is missing")
      if (allocated(error)) return
      call check(error, text == "unknown", "a state nothing assigned a spin to "// &
                 "should say so rather than be given one")

      call json%destroy()
      call data%destroy()
   end subroutine test_unrestricted_spin

end module test_mqc_json_writer

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_json_writer, only: collect_mqc_json_writer_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0
   testsuites = [new_testsuite("mqc_json_writer", collect_mqc_json_writer_tests)]

   do is = 1, size(testsuites)
      write (error_unit, fmt) "Testing:", testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, error_unit, stat)
   end do

   if (stat > 0) then
      write (error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
   end if
end program tester
