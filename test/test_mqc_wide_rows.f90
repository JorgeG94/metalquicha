module test_mqc_wide_rows
   !! Term-list rows wider than the MBE level
   !!
   !! A counterpoise row names its ghosts as well as its real monomers, and a
   !! full-cluster-basis row ghosts the whole system, so it can be wider than
   !! `max_level` and wider than `MAX_MBE_LEVEL`. Every consumer of the rows used
   !! to assume otherwise: a key held in a stack array of `MAX_MBE_LEVEL`, a copy
   !! of the first `max_level` columns, a table that listed rows by their total
   !! entries. These cases run such a list through each of them.
   !!
   !! The lists are the ones the scheme defines, built by hand in `wide_rows`
   !! from the definition: one own-basis row per monomer, and at level 2 one
   !! full-basis row per monomer and per pair, the rest of the system ghosted.
   !! Twelve monomers is past the stack bound of ten, which is the point. The
   !! energies are made up and different for every ghost pattern, so a term
   !! looked up in the wrong basis moves the total.
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use pic_types, only: dp, int64
   use pic_logger, only: logger => global_logger, error_level, verbose_level
   use mqc_mbe, only: compute_mbe
   use mqc_mbe_io, only: print_detailed_breakdown
   use mqc_result_types, only: calculation_result_t, mbe_result_t, SCF_NOT_CONVERGED
   use mqc_json_output_types, only: json_output_data_t
   use mqc_json_writer, only: write_json_output
   use mqc_fragment_table_writer, only: fragment_table_filename
   use mqc_io_helpers, only: set_output_json_filename, get_output_json_filename
   use mqc_frag_utils, only: generate_mbe_term_list
   use mqc_config_adapter, only: driver_config_t
   use mqc_combinatorics, only: COUNTERPOISE_VMFC, COUNTERPOISE_SSFC
   use mqc_physical_fragment, only: system_geometry_t
   use mqc_program_limits, only: MAX_MBE_LEVEL
   use json_module, only: json_file
   implicit none
   private
   public :: collect_mqc_wide_rows_tests

   real(dp), parameter :: TOLERANCE = 1.0e-12_dp
   real(dp), parameter :: RELATIVE_TOLERANCE = 1.0e-11_dp
      !! Of the total, for the hundred-monomer sum: its terms are large and the
      !! recursion and the check add them in a different order
   integer, parameter :: N_WIDE = 12
      !! Monomers in the wide fixtures; past `MAX_MBE_LEVEL`
   character(len=*), parameter :: LOG_PATH = "wide_rows_log.tmp"

contains

   subroutine collect_mqc_wide_rows_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("a_vmfc_table_padded_past_the_stack_bound_sums_unchanged", test_padded_vmfc), &
                  new_unittest("ghosted_rows_wider_than_the_stack_bound_are_subtracted", test_wide_expansion), &
                  new_unittest("the_written_tables_keep_the_ghosts_and_the_width", test_written_tables), &
                  new_unittest("the_detailed_breakdown_lists_wide_rows_at_their_real_level", &
                               test_detailed_breakdown) &
                  ]
   end subroutine collect_mqc_wide_rows_tests

   !---------------------------------------------------------------------------
   ! Fixtures
   !---------------------------------------------------------------------------

   subroutine make_chain(sys_geom, n_monomers)
      !! `n_monomers` one-atom monomers in a line, far enough apart for no screening
      type(system_geometry_t), intent(out) :: sys_geom
      integer, intent(in) :: n_monomers

      integer :: i

      sys_geom%n_monomers = n_monomers
      sys_geom%atoms_per_monomer = 1
      sys_geom%total_atoms = n_monomers
      sys_geom%charge = 0
      sys_geom%multiplicity = 1
      allocate (sys_geom%element_numbers(n_monomers))
      allocate (sys_geom%coordinates(3, n_monomers))
      sys_geom%element_numbers = 10  ! neon: closed shell, never bonded
      sys_geom%coordinates = 0.0_dp
      do i = 1, n_monomers
         sys_geom%coordinates(1, i) = real(i - 1, dp)*4.0_dp
      end do
   end subroutine make_chain

   function full_basis_row(members, n_monomers) result(row)
      !! The row of `members` real and every other monomer of the system ghosted
      integer, intent(in) :: members(:)
      integer, intent(in) :: n_monomers
      integer :: row(n_monomers)

      integer :: m, k

      row = 0
      row(1:size(members)) = members
      k = size(members)
      do m = 1, n_monomers
         if (any(members == m)) cycle
         k = k + 1
         row(k) = -m
      end do
   end function full_basis_row

   subroutine wide_rows(n_monomers, polymers, count)
      !! The full-cluster-basis list at level 2, `n_monomers` wide
      integer, intent(in) :: n_monomers
      integer, allocatable, intent(out) :: polymers(:, :)
      integer(int64), intent(out) :: count

      integer :: i, j

      allocate (polymers(2*n_monomers + n_monomers*(n_monomers - 1)/2, n_monomers))
      polymers = 0
      count = 0_int64
      do i = 1, n_monomers
         count = count + 1_int64
         polymers(count, 1) = i
      end do
      do i = 1, n_monomers
         count = count + 1_int64
         polymers(count, :) = full_basis_row([i], n_monomers)
      end do
      do i = 1, n_monomers
         do j = i + 1, n_monomers
            count = count + 1_int64
            polymers(count, :) = full_basis_row([i, j], n_monomers)
         end do
      end do
   end subroutine wide_rows

   function made_up(row) result(e)
      !! An energy for a row, symmetric in it and different for every ghost pattern
      integer, intent(in) :: row(:)
      real(dp) :: e

      integer :: a, b

      e = 0.0_dp
      do a = 1, size(row)
         if (row(a) == 0) cycle
         if (row(a) > 0) then
            e = e - 10.0_dp - 0.1_dp*real(row(a), dp)
         else
            e = e - 1.3_dp + 0.07_dp*real(row(a), dp)
         end if
         do b = a + 1, size(row)
            if (row(b) == 0) cycle
            e = e + 0.01_dp*real(mod(abs(row(a)*row(b)), 7) + 1, dp) &
                *(1.0_dp + 0.5_dp*real(count([row(a), row(b)] < 0), dp))
         end do
      end do
   end function made_up

   subroutine fill_results(polymers, count, results)
      !! Made-up energies for every row, and dipoles that are linear in them
      integer, intent(in) :: polymers(:, :)
      integer(int64), intent(in) :: count
      type(calculation_result_t), allocatable, intent(out) :: results(:)

      integer(int64) :: i

      allocate (results(count))
      do i = 1_int64, count
         results(i)%has_energy = .true.
         results(i)%energy%scf = made_up(polymers(i, :))
         allocate (results(i)%dipole(3))
         results(i)%dipole = [1.0_dp, 2.0_dp, -0.5_dp]*results(i)%energy%scf
         results(i)%has_dipole = .true.
      end do
   end subroutine fill_results

   !---------------------------------------------------------------------------
   ! The recursion
   !---------------------------------------------------------------------------

   subroutine test_padded_vmfc(error)
      !! A VMFC table padded with zero columns to twelve sums as the narrow one does
      !!
      !! The padding changes the width of every row and nothing else, so any
      !! consumer that read the width as the level, or held a row in a fixed
      !! buffer, would move the total.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      type(driver_config_t) :: config
      integer, allocatable :: polymers(:, :), padded(:, :)
      integer(int64) :: count
      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: narrow, wide

      call make_chain(sys_geom, N_WIDE)
      config%nlevel = 2
      config%counterpoise = "vmfc"
      call generate_mbe_term_list(sys_geom, config, 2, polymers, count)
      call check(error, size(polymers, 2) == 2, "the VMFC list should be level-2 wide")
      if (allocated(error)) return

      allocate (padded(size(polymers, 1), N_WIDE))
      padded = 0
      padded(:, 1:2) = polymers
      call check(error, size(padded, 2) > MAX_MBE_LEVEL, "the padded table must be wider than the stack bound")
      if (allocated(error)) return

      call fill_results(polymers, count, results)
      call narrow%allocate_dipole()
      call compute_mbe(polymers, count, 2, results, narrow, sys_geom=sys_geom, &
                       counterpoise_scheme=COUNTERPOISE_VMFC)
      call wide%allocate_dipole()
      call compute_mbe(padded, count, 2, results, wide, sys_geom=sys_geom, &
                       counterpoise_scheme=COUNTERPOISE_VMFC)

      call check(error, wide%total_energy, narrow%total_energy, thr=TOLERANCE, &
                 message="zero padding moved the VMFC total")
      if (allocated(error)) return
      call check(error, maxval(abs(wide%dipole - narrow%dipole)) < TOLERANCE, &
                 "zero padding moved the VMFC dipole")
      call narrow%destroy()
      call wide%destroy()
   end subroutine test_padded_vmfc

   subroutine test_wide_expansion(error)
      !! Rows that ghost most of the system are subtracted from their parents
      !!
      !! The subset key of a row names every ghost it carries, so it is as long
      !! as the system. Held in an array of `MAX_MBE_LEVEL` it ran off the end:
      !! unnoticed at twelve monomers, where a neighbouring local is overwritten
      !! and the answer survives, and fatally at a hundred. Twelve is the case
      !! that matters in practice and a hundred is the one that cannot pass by
      !! luck. The expected total is inclusion-exclusion done directly over
      !! full-basis energies and shares nothing with the recursion.
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: SIZES(2) = [N_WIDE, 100]
      integer :: k

      do k = 1, size(SIZES)
         call expand_and_compare(error, SIZES(k))
         if (allocated(error)) return
      end do
   end subroutine test_wide_expansion

   subroutine expand_and_compare(error, n_monomers)
      !! The level-2 full-basis expansion of `n_monomers` monomers against the direct sum
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: n_monomers

      type(system_geometry_t) :: sys_geom
      integer, allocatable :: polymers(:, :)
      integer(int64) :: count
      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: mbe_result
      real(dp) :: expected
      integer :: i, j
      character(len=16) :: label

      write (label, "(i0,a)") n_monomers, " monomers:"
      call make_chain(sys_geom, n_monomers)
      call wide_rows(n_monomers, polymers, count)
      call check(error, size(polymers, 2) > MAX_MBE_LEVEL, "the fixture must be wider than the stack bound")
      if (allocated(error)) return

      expected = 0.0_dp
      do i = 1, n_monomers
         expected = expected + made_up([i])
      end do
      do i = 1, n_monomers
         do j = i + 1, n_monomers
            expected = expected + made_up(full_basis_row([i, j], n_monomers)) &
                       - made_up(full_basis_row([i], n_monomers)) &
                       - made_up(full_basis_row([j], n_monomers))
         end do
      end do

      call fill_results(polymers, count, results)
      call mbe_result%allocate_dipole()
      call compute_mbe(polymers, count, 2, results, mbe_result, sys_geom=sys_geom, &
                       counterpoise_scheme=COUNTERPOISE_SSFC)

      call check(error, mbe_result%total_energy, expected, thr=RELATIVE_TOLERANCE*max(1.0_dp, abs(expected)), &
                 message=trim(label)//" the expansion is not full-basis inclusion-exclusion")
      if (allocated(error)) return
      call check(error, maxval(abs(mbe_result%dipole - [1.0_dp, 2.0_dp, -0.5_dp]*expected)) < &
                 2.0_dp*RELATIVE_TOLERANCE*max(1.0_dp, abs(expected)), &
                 trim(label)//" the dipole, linear in the energy, did not follow it")
      call mbe_result%destroy()
   end subroutine expand_and_compare

   !---------------------------------------------------------------------------
   ! The written tables
   !---------------------------------------------------------------------------

   subroutine test_written_tables(error)
      !! The JSON and the CSV carry the whole row, ghosts included, at full width
      !!
      !! One of the rows is marked unconverged, so the list of failures is
      !! checked for the row it should name rather than the first `max_level`
      !! columns of it.
      type(error_type), allocatable, intent(out) :: error

      type(system_geometry_t) :: sys_geom
      integer, allocatable :: polymers(:, :)
      integer(int64) :: count
      integer(int64), parameter :: FAILED_ROW = 30_int64
      type(calculation_result_t), allocatable :: results(:)
      type(mbe_result_t) :: mbe_result
      type(json_output_data_t) :: data
      type(json_file) :: json
      integer, allocatable :: indices(:), ghosts(:), named(:)
      character(len=2048) :: line, header
      integer :: unit, ios, i
      logical :: found, seen_wide_row

      call make_chain(sys_geom, N_WIDE)
      call wide_rows(N_WIDE, polymers, count)
      call fill_results(polymers, count, results)
      results(FAILED_ROW)%scf_status = SCF_NOT_CONVERGED
      call mbe_result%allocate_dipole()
      call compute_mbe(polymers, count, 2, results, mbe_result, sys_geom=sys_geom, &
                       json_data=data, counterpoise_scheme=COUNTERPOISE_SSFC)

      ! What compute_mbe hands the writers
      call check(error, size(data%polymers, 2) == N_WIDE, &
                 "the rows handed to the writers were cut to the level")
      if (allocated(error)) return
      call check(error, all(data%polymers == polymers(1:count, :)), "a row changed on its way to the writers")
      if (allocated(error)) return
      call check(error, allocated(data%unconverged_monomers), "the failed row was not collected")
      if (allocated(error)) return
      call check(error, size(data%unconverged_monomers, 2) == N_WIDE, &
                 "the failure list holds rows narrower than the list")
      if (allocated(error)) return
      call check(error, all(data%unconverged_monomers(1, :) == polymers(FAILED_ROW, :)), &
                 "the failure list lost the ghosts of the row that failed")
      if (allocated(error)) return

      ! The CSV, one column per entry of the row
      call set_output_json_filename("jw_wide.json")
      data%fragment_breakdown = "csv"
      call write_json_output(data)
      open (newunit=unit, file=trim(fragment_table_filename()), status="old", action="read", iostat=ios)
      call check(error, ios == 0, "the fragment table was not written")
      if (allocated(error)) return
      read (unit, "(a)") header
      call check(error, index(header, ",m12,") > 0, "the table has no column for the twelfth entry")
      if (allocated(error)) return
      call check(error, index(header, ",m13") == 0, "the table has a column past the row width")
      if (allocated(error)) return
      seen_wide_row = .false.
      do
         read (unit, "(a)", iostat=ios) line
         if (ios /= 0) exit
         call check(error, count_commas(line) == count_commas(header), &
                    "a row of the table is not as wide as its header")
         if (allocated(error)) exit
         if (index(line, ",2,1,2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,") > 0) seen_wide_row = .true.
      end do
      close (unit)
      if (allocated(error)) return
      call check(error, seen_wide_row, "the pair [1,2] with its ten ghosts is not in the table as it ran")
      if (allocated(error)) return

      ! The JSON: real monomers in `indices`, ghosts beside them
      data%fragment_breakdown = "json"
      call write_json_output(data)
      call json%initialize()
      call json%load_file(trim(get_output_json_filename()))

      call json%get("jw_wide.levels(2).fragments(1).indices", indices, found)
      call check(error, found, "a dimer has no indices")
      if (allocated(error)) return
      call json%get("jw_wide.levels(2).fragments(1).ghosts", ghosts, found)
      call check(error, found, "a full-basis dimer has no ghosts")
      if (allocated(error)) return
      call check(error, size(indices) == 2 .and. size(ghosts) == N_WIDE - 2, &
                 "a dimer should name two real monomers and ten ghosts")
      if (allocated(error)) return
      named = [indices, ghosts]
      call sort_ascending(named)
      call check(error, all(named == [(i, i=1, N_WIDE)]), "the real and ghosted monomers are not the whole system")
      if (allocated(error)) return

      ! Row 13 of the monomer level is the first full-basis monomer, [1,-2..-12]
      call json%get("jw_wide.levels(1).fragments(13).ghosts", ghosts, found)
      call check(error, found, "a full-basis monomer has no ghosts")
      if (allocated(error)) return
      call check(error, size(ghosts) == N_WIDE - 1, "a full-basis monomer should ghost the other eleven")
      if (allocated(error)) return

      ! An own-basis monomer is written exactly as it always was
      call json%info("jw_wide.levels(1).fragments(1).ghosts", found=found)
      call check(error,.not. found, "an own-basis monomer was given a ghosts key")
      if (allocated(error)) return

      call json%destroy()
      call data%destroy()
      call mbe_result%destroy()
      call remove_file(trim(get_output_json_filename()))
      call remove_file(trim(fragment_table_filename()))
   end subroutine test_written_tables

   subroutine remove_file(path)
      !! Delete a file this test wrote
      character(len=*), intent(in) :: path

      integer :: unit, ios

      open (newunit=unit, file=path, status="old", action="readwrite", iostat=ios)
      if (ios == 0) close (unit, status="delete")
   end subroutine remove_file

   function count_commas(text) result(n)
      character(len=*), intent(in) :: text
      integer :: n

      integer :: i

      n = 0
      do i = 1, len_trim(text)
         if (text(i:i) == ",") n = n + 1
      end do
   end function count_commas

   subroutine sort_ascending(values)
      !! Insertion sort of a handful of integers
      integer, intent(inout) :: values(:)

      integer :: i, j, held

      do i = 2, size(values)
         held = values(i)
         j = i - 1
         do while (j >= 1)
            if (values(j) <= held) exit
            values(j + 1) = values(j)
            j = j - 1
         end do
         values(j + 1) = held
      end do
   end subroutine sort_ascending

   !---------------------------------------------------------------------------
   ! The log
   !---------------------------------------------------------------------------

   subroutine test_detailed_breakdown(error)
      !! The verbose breakdown lists every row, at the level of its real monomers
      !!
      !! Every row of the fixture names twelve monomers, real or ghosted. Grouped
      !! by that total the table found no row at level 1 or 2 and listed nothing;
      !! grouped by the real monomers it has 24 monomers and 66 dimers.
      type(error_type), allocatable, intent(out) :: error

      integer, allocatable :: polymers(:, :)
      integer(int64) :: count
      real(dp), allocatable :: energies(:), deltas(:)
      character(len=2048) :: line
      integer :: unit, ios, console_level, n_rows_listed
      logical :: saw_monomers, saw_dimers, saw_wide_row

      call wide_rows(N_WIDE, polymers, count)
      allocate (energies(count), deltas(count))
      energies = -1.0_dp
      deltas = 0.0_dp

      ! Console off, file at verbose, so the lines can be read back
      call logger%configuration(level=console_level)
      call logger%configure(level=error_level)
      call logger%configure_file_output(LOG_PATH, level=verbose_level)
      call print_detailed_breakdown(polymers, count, 2, energies, deltas)
      call logger%close_log_file()
      call logger%configure(level=console_level)

      open (newunit=unit, file=LOG_PATH, status="old", action="read", iostat=ios)
      call check(error, ios == 0, "the log was not written")
      if (allocated(error)) return
      saw_monomers = .false.
      saw_dimers = .false.
      saw_wide_row = .false.
      n_rows_listed = 0
      do
         read (unit, "(a)", iostat=ios) line
         if (ios /= 0) exit
         if (index(line, "Monomers (24 fragments):") > 0) saw_monomers = .true.
         if (index(line, "Dimers (66 fragments):") > 0) saw_dimers = .true.
         if (index(line, "Fragment [") > 0) n_rows_listed = n_rows_listed + 1
         if (index(line, "Fragment [1,2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12]") > 0) saw_wide_row = .true.
      end do
      close (unit, status="delete")

      call check(error, saw_monomers, "the 24 monomer-level rows were not listed as monomers")
      if (allocated(error)) return
      call check(error, saw_dimers, "the 66 dimer-level rows were not listed as dimers")
      if (allocated(error)) return
      call check(error, n_rows_listed == 90, "not every row of the table was listed")
      if (allocated(error)) return
      call check(error, saw_wide_row, "a row was not printed with all twelve entries")
   end subroutine test_detailed_breakdown

end module test_mqc_wide_rows

program tester
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite, new_testsuite, testsuite_type
   use test_mqc_wide_rows, only: collect_mqc_wide_rows_tests
   implicit none
   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
                new_testsuite("mqc_wide_rows", collect_mqc_wide_rows_tests) &
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
