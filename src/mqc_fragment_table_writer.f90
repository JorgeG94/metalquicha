!! Flat per-fragment table output for MBE runs
module mqc_fragment_table_writer
   !! Writes the per-fragment MBE breakdown as a CSV sidecar next to the summary JSON.
   !!
   !! The breakdown is a homogeneous table of numbers, one row per fragment, which is
   !! the shape JSON handles worst: a repeated key string for every field of every
   !! row. Embedding it cost 18.3 us per fragment -- about seven minutes on a
   !! twenty-million-fragment run, and by far the largest part of the serial tail on
   !! rank 0 -- against 2.2 us for the same data written flat.
   !!
   !! It is also the shape the consumer wants. Convergence analysis and cross-method
   !! comparison mean pandas, and json.load has to hold the whole document resident:
   !! roughly 16 GB for twenty million fragments, times one per method being compared.
   !! A CSV streams, chunks, and joins.
   !!
   !! The summary -- total energy, per-level sums, thermochemistry, norms -- stays in
   !! the JSON, where a human can still read it.
   use pic_types, only: int64, dp
   use pic_logger, only: logger => global_logger
   use pic_timer, only: timer_type
   use pic_io, only: to_char
   use mqc_json_output_types, only: json_output_data_t
   use mqc_result_types, only: scf_status_label
   use mqc_physical_constants, only: HARTREE_TO_EV
   use mqc_io_helpers, only: get_basename
   implicit none
   private

   public :: write_fragment_table  !! Write per-fragment breakdown as CSV
   public :: fragment_table_filename  !! Sidecar name for a given run

contains

   function fragment_table_filename() result(filename)
      !! Sidecar name derived from the JSON output name: output_w1.json -> output_w1_fragments.csv
      character(len=256) :: filename
      filename = "output_"//trim(get_basename())//"_fragments.csv"
   end function fragment_table_filename

   subroutine write_fragment_table(data)
      !! Write one row per fragment: identity, energy, many-body correction, distance
      !!
      !! Monomer indices go out as fixed columns m1..m<max_level>, zero-filled, rather
      !! than a packed list. That keeps the file rectangular, and (level, m1..mL) is a
      !! stable key for joining the same system computed with different methods.
      type(json_output_data_t), intent(in) :: data

      integer :: unit, ios, j, level
      integer(int64) :: i
      logical :: have_scf
      logical :: have_orbitals
      logical :: have_energy, have_delta, have_distance
      type(timer_type) :: table_timer
      character(len=256) :: filename
      character(len=32) :: col
      character(len=64) :: row_fmt

      if (.not. allocated(data%polymers)) return
      if (data%fragment_count <= 0_int64) return

      filename = fragment_table_filename()
      open (newunit=unit, file=trim(filename), status="replace", action="write", iostat=ios)
      if (ios /= 0) then
         call logger%error("Could not open fragment table for writing: "//trim(filename))
         return
      end if

      call table_timer%start()

      ! Header
      write (unit, "(a)", advance="no") "frag_index,level"
      do j = 1, data%max_level
         write (col, "(a,i0)") ",m", j
         write (unit, "(a)", advance="no") trim(col)
      end do
      write (unit, "(a)") ",energy,delta_energy,distance,scf,homo,lumo,gap_ev"

      ! Presence of the value columns is fixed for the whole run, so decide once
      ! rather than per row.
      have_energy = allocated(data%fragment_energies)
      have_delta = allocated(data%delta_energies)
      have_distance = allocated(data%fragment_distances)
      have_scf = allocated(data%fragment_scf_status)
      have_orbitals = allocated(data%fragment_homo) .and. allocated(data%fragment_lumo) &
                      .and. allocated(data%fragment_has_orbitals)
      have_orbitals = allocated(data%fragment_homo) .and. allocated(data%fragment_lumo)

      ! Explicit repeat count for the monomer columns rather than an unlimited "*"
      ! group: the unlimited form emits the separator before it discovers the data is
      ! exhausted, which leaves a trailing comma and an extra column on every row.
      write (row_fmt, "(a,i0,a)") '(i0,",",i0,', data%max_level, '(",",i0))'

      ! Two write statements per row rather than one per column: statement overhead
      ! dominates at this row count, and collapsing them more than halves the cost.
      ! Reals go out at full precision, not a rounded display value -- screening
      ! studies difference the distances against the cutoffs that produced the list.
      do i = 1_int64, data%fragment_count
         level = count(data%polymers(i, :) > 0)

         write (unit, trim(row_fmt), advance="no") &
            i, level, (data%polymers(i, j), j=1, data%max_level)

         if (have_energy .and. have_delta .and. have_distance) then
            write (unit, '(3(",",es24.16))', advance="no") &
               data%fragment_energies(i), data%delta_energies(i), data%fragment_distances(i)
         else
            call write_optional_value(unit, have_energy, data%fragment_energies, i, .false.)
            call write_optional_value(unit, have_delta, data%delta_energies, i, .false.)
            call write_optional_value(unit, have_distance, data%fragment_distances, i, .false.)
         end if
         ! A word rather than the integer: this column is read by eye as often
         ! as by program, and "NO" in a spreadsheet is harder to scroll past
         ! than a 2.
         if (have_scf) then
            write (unit, "(a)", advance="no") ","//trim(scf_status_label(data%fragment_scf_status(i)))
         else
            write (unit, "(a)", advance="no") ",?"
         end if

         ! Repeated in eV because that is the unit a gap gets compared in, and
         ! converting a column of Hartrees by hand is how a factor of 27 ends
         ! up in a plot. Blank, not zero, where the method reported no pair:
         ! a gap of zero is a claim about the fragment.
         if (have_orbitals .and. data%fragment_has_orbitals(i)) then
            write (unit, '(2(",",es24.16),",",f12.6)') data%fragment_homo(i), &
               data%fragment_lumo(i), &
               (data%fragment_lumo(i) - data%fragment_homo(i))*HARTREE_TO_EV
         else
            write (unit, "(a)") ",,,"
         end if
      end do

      close (unit)
      call table_timer%stop()

      call logger%info("Fragment table written to "//trim(filename)//" ("// &
                       to_char(data%fragment_count)//" rows in "// &
                       to_char(table_timer%get_elapsed_time())//" s)")

   end subroutine write_fragment_table

   subroutine write_optional_value(unit, present_col, values, i, last)
      !! Emit one comma-separated value, or an empty cell if the column is absent
      integer, intent(in) :: unit
      logical, intent(in) :: present_col
      real(dp), allocatable, intent(in) :: values(:)
      integer(int64), intent(in) :: i
      logical, intent(in) :: last  !! .true. ends the record

      if (present_col) then
         if (last) then
            write (unit, '(",",es24.16)') values(i)
         else
            write (unit, '(",",es24.16)', advance="no") values(i)
         end if
      else
         if (last) then
            write (unit, "(a)") ","
         else
            write (unit, "(a)", advance="no") ","
         end if
      end if
   end subroutine write_optional_value

end module mqc_fragment_table_writer
