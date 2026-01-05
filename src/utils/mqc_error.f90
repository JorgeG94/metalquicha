!! Error handling module for metalquicha
!! Provides a unified error type to replace stat/errmsg pairs
!! Enhanced with stack trace support for better debugging
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

   !! Stack trace configuration
   integer, parameter :: MAX_STACK_DEPTH = 20
   integer, parameter :: MAX_LOCATION_LEN = 128

   !! Unified error type with stack trace support
   type :: error_t
      integer :: code = SUCCESS  !! Error code (0 = no error)
      character(len=:), allocatable :: message  !! Error message

      !! Stack trace support
      integer :: stack_depth = 0  !! Current stack depth
      character(len=MAX_LOCATION_LEN) :: call_stack(MAX_STACK_DEPTH)  !! Call locations
   contains
      procedure :: has_error => error_has_error
      procedure :: set => error_set
      procedure :: clear => error_clear
      procedure :: get_code => error_get_code
      procedure :: get_message => error_get_message
      procedure :: add_context => error_add_context
      procedure :: get_full_trace => error_get_full_trace
      procedure :: print_trace => error_print_trace
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
      !! Resets the stack trace
      class(error_t), intent(inout) :: this
      integer, intent(in) :: code
      character(len=*), intent(in) :: message

      this%code = code
      this%message = trim(message)
      this%stack_depth = 0  ! Reset stack when setting new error
   end subroutine error_set

   pure subroutine error_clear(this)
      !! Clear the error state and stack trace
      class(error_t), intent(inout) :: this
      this%code = SUCCESS
      this%stack_depth = 0
      if (allocated(this%message)) deallocate (this%message)
   end subroutine error_clear

   pure function error_get_code(this) result(code)
      !! Get the error code
      class(error_t), intent(in) :: this
      integer :: code
      code = this%code
   end function error_get_code

   pure function error_get_message(this) result(message)
      !! Get the error message (without stack trace)
      class(error_t), intent(in) :: this
      character(len=:), allocatable :: message
      if (allocated(this%message)) then
         message = this%message
      else
         message = ""
      end if
   end function error_get_message

   pure subroutine error_add_context(this, location)
      !! Add a call location to the stack trace
      !! Typically called when propagating errors upward
      !!
      !! Example:
      !!   call some_routine(..., error)
      !!   if (error%has_error()) then
      !!      call error%add_context("mqc_mbe:compute_energy")
      !!      return
      !!   end if
      class(error_t), intent(inout) :: this
      character(len=*), intent(in) :: location

      if (this%stack_depth < MAX_STACK_DEPTH) then
         this%stack_depth = this%stack_depth + 1
         this%call_stack(this%stack_depth) = location
      end if
      ! If stack is full, silently ignore (could print warning in non-pure version)
   end subroutine error_add_context

   function error_get_full_trace(this) result(trace)
      !! Get complete error message with stack trace
      !! Returns a multi-line string with error and call stack
      class(error_t), intent(in) :: this
      character(len=:), allocatable :: trace
      character(len=2048) :: buffer
      integer :: i, pos

      if (.not. this%has_error()) then
         trace = ""
         return
      end if

      ! Build error message
      write (buffer, '(A,I0,A)') "Error ", this%code, ": "
      pos = len_trim(buffer) + 1

      if (allocated(this%message)) then
         buffer(pos:) = this%message
         pos = len_trim(buffer) + 1
      end if

      ! Add stack trace if available
      if (this%stack_depth > 0) then
         buffer(pos:) = new_line('a')//"Call stack (most recent first):"
         pos = len_trim(buffer) + 1

         do i = this%stack_depth, 1, -1
            write (buffer(pos:), '(A,I0,A)') new_line('a')//"  [", i, "] "
            pos = len_trim(buffer) + 1
            buffer(pos:) = trim(this%call_stack(i))
            pos = len_trim(buffer) + 1
         end do
      end if

      trace = trim(buffer)
   end function error_get_full_trace

   subroutine error_print_trace(this, unit)
      !! Print error with stack trace to specified unit
      !! If unit not specified, prints to stdout (unit 6)
      class(error_t), intent(in) :: this
      integer, intent(in), optional :: unit
      integer :: out_unit, i

      out_unit = 6  ! stdout
      if (present(unit)) out_unit = unit

      if (.not. this%has_error()) return

      ! Print error message
      write (out_unit, '(A,I0,A)', advance='no') "Error ", this%code, ": "
      if (allocated(this%message)) then
         write (out_unit, '(A)') trim(this%message)
      else
         write (out_unit, '(A)') "(no message)"
      end if

      ! Print stack trace if available
      if (this%stack_depth > 0) then
         write (out_unit, '(A)') "Call stack (most recent first):"
         do i = this%stack_depth, 1, -1
            write (out_unit, '(A,I0,A)', advance='no') "  [", i, "] "
            write (out_unit, '(A)') trim(this%call_stack(i))
         end do
      end if
   end subroutine error_print_trace

end module mqc_error
