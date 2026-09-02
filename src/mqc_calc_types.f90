!! Calculation type constants for quantum chemistry calculations
module mqc_calc_types
   !! The `CALC_TYPE_*` constants, and conversion between them and the driver
   !! names a deck spells.
   use pic_types, only: int32
   implicit none
   private

   ! Public constants
   public :: CALC_TYPE_ENERGY, CALC_TYPE_GRADIENT, CALC_TYPE_HESSIAN
   public :: CALC_TYPE_MAKEFP
   public :: CALC_TYPE_OPTIMIZE
   public :: CALC_TYPE_CONFORMERS
   public :: CALC_TYPE_UNKNOWN

   ! Public functions
   public :: calc_type_from_string, calc_type_to_string

   ! Calculation type constants
   integer(int32), parameter :: CALC_TYPE_UNKNOWN = 0
   integer(int32), parameter :: CALC_TYPE_ENERGY = 1
   integer(int32), parameter :: CALC_TYPE_GRADIENT = 2
   integer(int32), parameter :: CALC_TYPE_HESSIAN = 3
   integer(int32), parameter :: CALC_TYPE_MAKEFP = 4
   !! Build an effective fragment potential and write it, computing no energy.
   integer(int32), parameter :: CALC_TYPE_OPTIMIZE = 5
   !! Minimize the geometry, which is a loop over gradient calculations rather
   !! than a calculation. Handled above `run_calculation`, not inside it.
   integer(int32), parameter :: CALC_TYPE_CONFORMERS = 6
   !! Search for conformers, a sampling run wrapped around many gradients
   !! rather than a calculation. Driven above `run_calculation`, like OPTIMIZE,
   !! and not dispatchable from inside it.

contains

   pure function calc_type_from_string(calc_type_str) result(calc_type)
      !! Calculation type for a driver name
      !!
      !! Case-insensitive. `CALC_TYPE_UNKNOWN` for anything unrecognized.
      character(len=*), intent(in) :: calc_type_str  !! Input string (e.g., "energy", "gradient")
      integer(int32) :: calc_type                     !! Output integer constant

      character(len=len_trim(calc_type_str)) :: lower_str
      integer :: i

      ! Convert to lowercase for case-insensitive comparison
      lower_str = trim(adjustl(calc_type_str))
      do i = 1, len(lower_str)
         if (lower_str(i:i) >= "A" .and. lower_str(i:i) <= "Z") then
            lower_str(i:i) = achar(iachar(lower_str(i:i)) + 32)
         end if
      end do

      ! Match against known types
      select case (lower_str)
      case ("energy")
         calc_type = CALC_TYPE_ENERGY
      case ("gradient")
         calc_type = CALC_TYPE_GRADIENT
      case ("hessian")
         calc_type = CALC_TYPE_HESSIAN
      case ("makefp", "makeefp")
         ! GAMESS spells the run type MAKEFP; "makeefp" is the common variant.
         calc_type = CALC_TYPE_MAKEFP
      case ("optimize", "optimization", "opt")
         ! "optimize" is the QCSchema driver name the decks here use.
         calc_type = CALC_TYPE_OPTIMIZE
      case ("conformers", "conformer", "crest")
         calc_type = CALC_TYPE_CONFORMERS
      case default
         calc_type = CALC_TYPE_UNKNOWN
      end select

   end function calc_type_from_string

   pure function calc_type_to_string(calc_type) result(calc_type_str)
      !! The driver name a calculation type constant stands for
      integer(int32), intent(in) :: calc_type         !! Input integer constant
      character(len=:), allocatable :: calc_type_str  !! Output string representation

      select case (calc_type)
      case (CALC_TYPE_ENERGY)
         calc_type_str = "energy"
      case (CALC_TYPE_GRADIENT)
         calc_type_str = "gradient"
      case (CALC_TYPE_HESSIAN)
         calc_type_str = "hessian"
      case (CALC_TYPE_MAKEFP)
         calc_type_str = "makefp"
      case (CALC_TYPE_OPTIMIZE)
         calc_type_str = "optimize"
      case (CALC_TYPE_CONFORMERS)
         calc_type_str = "conformers"
      case default
         calc_type_str = "unknown"
      end select

   end function calc_type_to_string

end module mqc_calc_types
