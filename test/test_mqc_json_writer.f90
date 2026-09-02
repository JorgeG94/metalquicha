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
                  new_unittest("a_fingerprint_is_written_when_there_is_one", test_fingerprint) &
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
