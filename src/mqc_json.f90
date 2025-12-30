!! JSON output utilities for multi-molecule calculations
module mqc_json
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   implicit none
   private

   public :: merge_multi_molecule_json

contains

   subroutine merge_multi_molecule_json(individual_files, nmol)
      !! Merge individual molecule JSON files into a single combined file
      use mqc_io_helpers, only: get_molecule_name

      character(len=256), intent(in) :: individual_files(:)
      integer, intent(in) :: nmol

      integer :: imol, unit_in, unit_out, io_stat, slash_pos, dot_pos
      character(len=10000) :: line
      character(len=256) :: output_file, basename
      logical :: file_exists

      ! Determine combined output filename from first individual file
      ! Example: "output_multi_structure_molecule_1.json" -> "output_multi_structure.json"
      basename = individual_files(1)
      slash_pos = index(basename, '/', back=.true.)
      if (slash_pos > 0) then
         basename = basename(slash_pos + 1:)
      end if

      ! Remove "_molecule_1" or similar suffix
      dot_pos = index(basename, '_molecule_')
      if (dot_pos > 0) then
         output_file = basename(1:dot_pos - 1)//".json"
      else
         output_file = "output_combined.json"
      end if

      ! Open combined output file
      open (newunit=unit_out, file=trim(output_file), status='replace', action='write', iostat=io_stat)
      if (io_stat /= 0) then
         call logger%error("Failed to open "//trim(output_file)//" for writing")
         return
      end if

      call logger%info("Merging "//to_char(nmol)//" molecule JSON files into "//trim(output_file))

      ! Write opening brace and top-level key (basename without "output_" and ".json")
      dot_pos = index(output_file, '.json')
      if (dot_pos > 0) then
         basename = output_file(8:dot_pos - 1)  ! Skip "output_"
      else
         basename = "combined"
      end if

      write (unit_out, '(a)') "{"
      write (unit_out, '(a)') '  "'//trim(basename)//'": {'

      ! Process each individual JSON file
      do imol = 1, nmol
         inquire (file=trim(individual_files(imol)), exist=file_exists)
         if (.not. file_exists) cycle

         open (newunit=unit_in, file=trim(individual_files(imol)), status='old', action='read', iostat=io_stat)
         if (io_stat /= 0) cycle

         ! Read all lines from the individual JSON file
         call read_json_content(unit_in, imol, unit_out, individual_files(imol))

         close (unit_in)

         ! Delete individual file
         open (newunit=unit_in, file=trim(individual_files(imol)), status='old', action='readwrite')
         close (unit_in, status='delete')
      end do

      ! Close last molecule
      write (unit_out, '(a)') '    }'

      ! Close top-level key and file
      write (unit_out, '(a)') '  }'
      write (unit_out, '(a)') '}'

      close (unit_out)
      call logger%info("Combined JSON written to "//trim(output_file))

   end subroutine merge_multi_molecule_json

   subroutine read_json_content(unit_in, mol_index, unit_out, filename)
      !! Read and write JSON content from an individual molecule file
      !! Properly handles nested structures from fragmented calculations
      use mqc_io_helpers, only: get_molecule_name

      integer, intent(in) :: unit_in, mol_index, unit_out
      character(len=*), intent(in) :: filename

      character(len=10000), allocatable :: all_lines(:)
      character(len=10000) :: line
      integer :: io_stat, nlines, i

      ! Read all lines into memory
      allocate (all_lines(1000))  ! Reasonable size for most JSON files
      nlines = 0

      do
         read (unit_in, '(a)', iostat=io_stat) line
         if (io_stat /= 0) exit
         nlines = nlines + 1
         if (nlines > size(all_lines)) then
            ! Reallocate if needed
            call logger%error("JSON file too large: "//trim(filename))
            return
         end if
         all_lines(nlines) = line
      end do

      ! Lines structure:
      ! 1: "{"
      ! 2: '  "molecule_name": {'
      ! 3..(n-2): content
      ! n-1: "  }"
      ! n: "}"

      if (nlines < 3) then
         call logger%error("Invalid JSON structure: "//trim(filename))
         return
      end if

      ! Write molecule key (extracted from filename)
      if (mol_index > 1) write (unit_out, '(a)') '    },'
      write (unit_out, '(a)') '    "'//trim(get_molecule_name(filename))//'" : {'

      ! Write all content lines (from line 3 to line n-2)
      do i = 3, nlines - 2
         write (unit_out, '(a)') '  '//trim(all_lines(i))  ! Add 2 spaces for proper indentation
      end do

      deallocate (all_lines)

   end subroutine read_json_content

end module mqc_json
