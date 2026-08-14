!! Small string helpers that have to work inside `pure` procedures
module mqc_string_utils
   !! Formatting helpers for building messages and JSON paths.
   !!
   !! `pic_io`'s `to_char` already does this and is what most of the program
   !! uses, but it is not `pure`, so it cannot be called from the validators and
   !! readers that are. That is the whole reason six near-identical private
   !! copies of "integer to string" had accumulated -- as `int_to_string`,
   !! `int_text` and `integer_to_key` -- rather than anyone preferring their own.
   !!
   !! So: use `to_char` in ordinary code, and this where `pure` is required.
   use pic_types, only: dp
   implicit none
   private

   public :: int_to_text

contains

   pure function int_to_text(value) result(text)
      !! An integer as its shortest decimal string, no padding
      integer, intent(in) :: value
      character(len=:), allocatable :: text

      character(len=32) :: buffer

      write (buffer, "(I0)") value
      text = trim(buffer)
   end function int_to_text

end module mqc_string_utils
