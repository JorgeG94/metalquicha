!! IO helper utilities for file naming and string operations
!! Provides utilities for output filename management and string parsing
module mqc_io_helpers
   implicit none
   private

   character(len=256), save :: output_json_filename = "results.json"
   character(len=256), save :: current_basename = ""
   public :: set_output_json_filename, get_output_json_filename, get_basename
   public :: set_molecule_suffix
   public :: get_molecule_name, ends_with

contains

   subroutine set_output_json_filename(input_filename)
      !! Set the JSON output filename based on input filename
      !! Example: "water.mqc" -> "output_water.json"
      character(len=*), intent(in) :: input_filename
      integer :: dot_pos, slash_pos
      character(len=256) :: basename

      ! Find last slash (if any) to extract basename
      slash_pos = index(input_filename, '/', back=.true.)
      if (slash_pos > 0) then
         basename = input_filename(slash_pos + 1:)
      else
         basename = input_filename
      end if

      ! Find last dot to remove extension
      dot_pos = index(basename, '.', back=.true.)
      if (dot_pos > 0) then
         basename = basename(1:dot_pos - 1)
      end if

      ! Store basename for later use
      current_basename = trim(basename)

      ! Construct output filename: output_<basename>.json
      output_json_filename = "output_"//trim(basename)//".json"

   end subroutine set_output_json_filename

   subroutine set_molecule_suffix(suffix)
      !! Append a suffix to the output filename (e.g., for multi-molecule mode)
      !! Example: suffix="_mol1" -> "output_multi_structure_mol1.json"
      character(len=*), intent(in) :: suffix

      if (len_trim(current_basename) > 0) then
         output_json_filename = "output_"//trim(current_basename)//trim(suffix)//".json"
      end if

   end subroutine set_molecule_suffix

   function get_output_json_filename() result(filename)
      !! Get the current JSON output filename
      character(len=256) :: filename
      filename = trim(output_json_filename)
   end function get_output_json_filename

   function get_basename() result(basename)
      !! Get the base name without "output_" prefix and ".json" suffix
      !! Example: "output_w1.json" -> "w1"
      character(len=256) :: basename
      integer :: start_pos, end_pos

      ! Remove "output_" prefix (7 characters)
      start_pos = 8

      ! Find ".json" suffix
      end_pos = index(output_json_filename, '.json', back=.true.) - 1

      if (end_pos > start_pos) then
         basename = output_json_filename(start_pos:end_pos)
      else
         basename = "unknown"
      end if
   end function get_basename

   function get_molecule_name(filename) result(name)
      !! Extract molecule name from filename
      !! Example: "output_multi_structure_molecule_1.json" -> "molecule_1"
      character(len=*), intent(in) :: filename
      character(len=256) :: name
      integer :: start_pos, end_pos

      ! Find "_molecule_" or similar pattern
      start_pos = index(filename, '_molecule_')
      if (start_pos == 0) start_pos = index(filename, '_mol_')

      if (start_pos > 0) then
         start_pos = start_pos + 1  ! Skip leading underscore
         end_pos = index(filename, '.json') - 1
         if (end_pos > start_pos) then
            name = filename(start_pos:end_pos)
         else
            name = "unknown"
         end if
      else
         name = "unknown"
      end if
   end function get_molecule_name

   logical function ends_with(str, suffix)
      !! Check if string ends with suffix
      character(len=*), intent(in) :: str, suffix
      integer :: str_len, suffix_len

      str_len = len_trim(str)
      suffix_len = len_trim(suffix)

      if (suffix_len > str_len) then
         ends_with = .false.
         return
      end if

      ends_with = (str(str_len - suffix_len + 1:str_len) == suffix)
   end function ends_with

end module mqc_io_helpers
