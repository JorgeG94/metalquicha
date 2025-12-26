!! Calculation type constants for quantum chemistry calculations
module mqc_calc_types
   !! Defines integer constants for calculation types to avoid string comparisons
   !! throughout the codebase. Provides conversion utilities between string
   !! representations and integer constants.
   use pic_types, only: int32
   implicit none
   private

   ! Public constants
   public :: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   public :: CALC_TYPE_UNKNOWN

   ! Public functions
   public :: calc_type_from_string, calc_type_to_string

   ! Calculation type constants
   integer(int32), parameter :: CALC_TYPE_UNKNOWN = 0
   integer(int32), parameter :: CALC_TYPE_ENERGY = 1
   integer(int32), parameter :: CALC_TYPE_GRADIENT = 2
   integer(int32), parameter :: CALC_TYPE_HESSIAN = 3

contains

   pure function calc_type_from_string(calc_type_str) result(calc_type)
      !! Convert calculation type string to integer constant
      !!
      !! Performs case-insensitive comparison and returns appropriate constant.
      !! Returns CALC_TYPE_UNKNOWN for unrecognized strings.
      character(len=*), intent(in) :: calc_type_str  !! Input string (e.g., "energy", "gradient")
      integer(int32) :: calc_type                     !! Output integer constant

      character(len=len_trim(calc_type_str)) :: lower_str
      integer :: i

      ! Convert to lowercase for case-insensitive comparison
      lower_str = trim(adjustl(calc_type_str))
      do i = 1, len(lower_str)
         if (lower_str(i:i) >= 'A' .and. lower_str(i:i) <= 'Z') then
            lower_str(i:i) = achar(iachar(lower_str(i:i)) + 32)
         end if
      end do

      ! Match against known types
      select case (lower_str)
      case ('energy')
         calc_type = CALC_TYPE_ENERGY
      case ('gradient')
         calc_type = CALC_TYPE_GRADIENT
      case ('hessian')
         calc_type = CALC_TYPE_HESSIAN
      case default
         calc_type = CALC_TYPE_UNKNOWN
      end select

   end function calc_type_from_string

   pure function calc_type_to_string(calc_type) result(calc_type_str)
      !! Convert calculation type integer constant to string
      !!
      !! Provides human-readable string representation of calculation type.
      integer(int32), intent(in) :: calc_type         !! Input integer constant
      character(len=:), allocatable :: calc_type_str  !! Output string representation

      select case (calc_type)
      case (CALC_TYPE_ENERGY)
         calc_type_str = "energy"
      case (CALC_TYPE_GRADIENT)
         calc_type_str = "gradient"
      case (CALC_TYPE_HESSIAN)
         calc_type_str = "hessian"
      case default
         calc_type_str = "unknown"
      end select

   end function calc_type_to_string

end module mqc_calc_types
