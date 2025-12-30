!! Module to store the base name for output files
!! This is set once from the input filename and used by all JSON writers
module mqc_output_filename
   implicit none
   private

   character(len=256), save :: output_json_filename = "results.json"
   character(len=256), save :: current_basename = ""
   public :: set_output_json_filename, get_output_json_filename, get_basename
   public :: set_molecule_suffix

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

end module mqc_output_filename
