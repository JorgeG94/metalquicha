!! IO helper utilities for file naming and string operations
!! Provides utilities for output filename management and string parsing
module mqc_io_helpers
   use pic_ascii, only: to_lower
   implicit none
   private

   character(len=256) :: output_json_filename = "results.json"
   character(len=256) :: current_basename = ""

   public :: set_output_json_filename, get_output_json_filename, get_basename
   public :: set_molecule_suffix
   public :: get_molecule_name, ends_with

contains

   subroutine set_output_json_filename(input_filename)
      !! Set the JSON output filename based on input filename
      !! Example: "water.json" -> "output_water.json"
      character(len=*), intent(in) :: input_filename
      integer :: dot_pos, slash_pos
      character(len=256) :: basename

      ! Find last slash (if any) to extract basename
      slash_pos = index(input_filename, "/", back=.true.)
      if (slash_pos > 0) then
         basename = input_filename(slash_pos + 1:)
      else
         basename = input_filename
      end if

      ! Find last dot to remove extension
      dot_pos = index(basename, ".", back=.true.)
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
      end_pos = index(output_json_filename, ".json", back=.true.) - 1

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
      start_pos = index(filename, "_molecule_")
      if (start_pos == 0) start_pos = index(filename, "_mol_")

      if (start_pos > 0) then
         start_pos = start_pos + 1  ! Skip leading underscore
         end_pos = index(filename, ".json") - 1
         if (end_pos > start_pos) then
            name = filename(start_pos:end_pos)
         else
            name = "unknown"
         end if
      else
         name = "unknown"
      end if
   end function get_molecule_name

   pure function ends_with(text, suffix) result(matches)
      !! Case-insensitive suffix test, for choosing a format by filename
      !!
      !! Insensitive because every caller is matching a file extension, and an
      !! extension's case is the user's typing rather than a fact about the
      !! file: `WATER.JSON` is a JSON deck.
      character(len=*), intent(in) :: text, suffix
      logical :: matches

      integer :: n, m

      n = len_trim(text)
      m = len_trim(suffix)
      matches = .false.
      if (m > n .or. m < 1) return
      matches = (to_lower(text(n - m + 1:n)) == to_lower(suffix(1:m)))
   end function ends_with

end module mqc_io_helpers
