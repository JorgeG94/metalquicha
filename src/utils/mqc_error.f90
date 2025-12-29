!! Error handling module for metalquicha
!! Provides a unified error type to replace stat/errmsg pairs
module mqc_error
   implicit none
   private

   public :: error_t
   public :: SUCCESS, ERROR_GENERIC, ERROR_IO, ERROR_PARSE, ERROR_VALIDATION

   !! Error codes
   integer, parameter :: SUCCESS = 0
   integer, parameter :: ERROR_GENERIC = 1
   integer, parameter :: ERROR_IO = 2
   integer, parameter :: ERROR_PARSE = 3
   integer, parameter :: ERROR_VALIDATION = 4

   !! Unified error type
   type :: error_t
      integer :: code = SUCCESS  !! Error code (0 = no error)
      character(len=:), allocatable :: message  !! Error message
   contains
      procedure :: has_error => error_has_error
      procedure :: set => error_set
      procedure :: clear => error_clear
      procedure :: get_code => error_get_code
      procedure :: get_message => error_get_message
   end type error_t

contains

   pure function error_has_error(this) result(has_err)
      !! Check if an error is set
      class(error_t), intent(in) :: this
      logical :: has_err
      has_err = (this%code /= SUCCESS)
   end function error_has_error

   pure subroutine error_set(this, code, message)
      !! Set an error with code and message
      class(error_t), intent(inout) :: this
      integer, intent(in) :: code
      character(len=*), intent(in) :: message

      this%code = code
      this%message = trim(message)
   end subroutine error_set

   pure subroutine error_clear(this)
      !! Clear the error state
      class(error_t), intent(inout) :: this
      this%code = SUCCESS
      if (allocated(this%message)) deallocate (this%message)
   end subroutine error_clear

   pure function error_get_code(this) result(code)
      !! Get the error code
      class(error_t), intent(in) :: this
      integer :: code
      code = this%code
   end function error_get_code

   pure function error_get_message(this) result(message)
      !! Get the error message
      class(error_t), intent(in) :: this
      character(len=:), allocatable :: message
      if (allocated(this%message)) then
         message = this%message
      else
         message = ""
      end if
   end function error_get_message

end module mqc_error
