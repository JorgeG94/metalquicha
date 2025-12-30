!! Module to store the base name for output files
!! This is set once from the input filename and used by all JSON writers
module mqc_output_filename
   implicit none
   private

   character(len=256), save :: output_json_filename = "results.json"
   public :: set_output_json_filename, get_output_json_filename

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

      ! Construct output filename: output_<basename>.json
      output_json_filename = "output_"//trim(basename)//".json"

   end subroutine set_output_json_filename

   function get_output_json_filename() result(filename)
      !! Get the current JSON output filename
      character(len=256) :: filename
      filename = trim(output_json_filename)
   end function get_output_json_filename

end module mqc_output_filename
