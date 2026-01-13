!! Method type constants for quantum chemistry methods
module mqc_method_types
   !! Defines integer constants for quantum chemistry methods to avoid string comparisons
   !! throughout the codebase. Provides conversion utilities between string
   !! representations and integer constants.
   use pic_types, only: int32
   implicit none
   private

   ! Public constants
   public :: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2, METHOD_TYPE_HF
   public :: METHOD_TYPE_DFT, METHOD_TYPE_MCSCF
   public :: METHOD_TYPE_UNKNOWN

   ! Public functions
   public :: method_type_from_string, method_type_to_string

   ! Method type constants
   integer(int32), parameter :: METHOD_TYPE_UNKNOWN = 0
   integer(int32), parameter :: METHOD_TYPE_GFN1 = 1
   integer(int32), parameter :: METHOD_TYPE_GFN2 = 2
   integer(int32), parameter :: METHOD_TYPE_HF = 3
   integer(int32), parameter :: METHOD_TYPE_DFT = 4
   integer(int32), parameter :: METHOD_TYPE_MCSCF = 5

contains

   pure function method_type_from_string(method_str) result(method_type)
      !! Convert method type string to integer constant
      !!
      !! Performs case-insensitive comparison and returns appropriate constant.
      !! Returns METHOD_TYPE_UNKNOWN for unrecognized strings.
      character(len=*), intent(in) :: method_str  !! Input string (e.g., "gfn1", "gfn2", "hf")
      integer(int32) :: method_type                !! Output integer constant

      character(len=len_trim(method_str)) :: lower_str
      integer :: i

      ! Convert to lowercase for case-insensitive comparison
      lower_str = trim(adjustl(method_str))
      do i = 1, len(lower_str)
         if (lower_str(i:i) >= 'A' .and. lower_str(i:i) <= 'Z') then
            lower_str(i:i) = achar(iachar(lower_str(i:i)) + 32)
         end if
      end do

      ! Match against known types
      select case (lower_str)
      case ('gfn1')
         method_type = METHOD_TYPE_GFN1
      case ('gfn2')
         method_type = METHOD_TYPE_GFN2
      case ('hf', 'rhf', 'uhf')
         method_type = METHOD_TYPE_HF
      case ('dft', 'ks', 'kohn-sham')
         method_type = METHOD_TYPE_DFT
      case ('mcscf', 'casscf', 'casci')
         method_type = METHOD_TYPE_MCSCF
      case default
         method_type = METHOD_TYPE_UNKNOWN
      end select

   end function method_type_from_string

   pure function method_type_to_string(method_type) result(method_str)
      !! Convert method type integer constant to string
      !!
      !! Provides human-readable string representation of method type.
      integer(int32), intent(in) :: method_type         !! Input integer constant
      character(len=:), allocatable :: method_str       !! Output string representation

      select case (method_type)
      case (METHOD_TYPE_GFN1)
         method_str = "gfn1"
      case (METHOD_TYPE_GFN2)
         method_str = "gfn2"
      case (METHOD_TYPE_HF)
         method_str = "hf"
      case (METHOD_TYPE_DFT)
         method_str = "dft"
      case (METHOD_TYPE_MCSCF)
         method_str = "mcscf"
      case default
         method_str = "unknown"
      end select

   end function method_type_to_string

end module mqc_method_types
