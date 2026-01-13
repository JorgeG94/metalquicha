!! Method type constants for quantum chemistry methods
module mqc_method_types
   !! Defines integer constants for quantum chemistry methods to avoid string comparisons
   !! throughout the codebase. Provides conversion utilities between string
   !! representations and integer constants.
   use pic_types, only: int32
   implicit none
   private

   ! Public constants - Semi-empirical
   public :: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
   ! Public constants - SCF methods
   public :: METHOD_TYPE_HF, METHOD_TYPE_DFT
   ! Public constants - Multi-reference
   public :: METHOD_TYPE_MCSCF
   ! Public constants - Correlation methods
   public :: METHOD_TYPE_MP2, METHOD_TYPE_CCSD, METHOD_TYPE_CCSD_T
   public :: METHOD_TYPE_MP2_F12, METHOD_TYPE_CCSD_F12, METHOD_TYPE_CCSD_T_F12
   public :: METHOD_TYPE_UNKNOWN

   ! Public functions
   public :: method_type_from_string, method_type_to_string

   ! Method type constants
   integer(int32), parameter :: METHOD_TYPE_UNKNOWN = 0

   ! Semi-empirical (1-9)
   integer(int32), parameter :: METHOD_TYPE_GFN1 = 1
   integer(int32), parameter :: METHOD_TYPE_GFN2 = 2

   ! SCF methods (10-19)
   integer(int32), parameter :: METHOD_TYPE_HF = 10
   integer(int32), parameter :: METHOD_TYPE_DFT = 11

   ! Multi-reference (20-29)
   integer(int32), parameter :: METHOD_TYPE_MCSCF = 20

   ! Perturbation theory (30-39)
   integer(int32), parameter :: METHOD_TYPE_MP2 = 30
   integer(int32), parameter :: METHOD_TYPE_MP2_F12 = 31

   ! Coupled cluster (40-59)
   integer(int32), parameter :: METHOD_TYPE_CCSD = 40
   integer(int32), parameter :: METHOD_TYPE_CCSD_T = 41      !! CCSD(T)
   integer(int32), parameter :: METHOD_TYPE_CCSD_F12 = 42
   integer(int32), parameter :: METHOD_TYPE_CCSD_T_F12 = 43  !! CCSD(T)-F12

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
         ! Semi-empirical
      case ('gfn1', 'gfn1-xtb')
         method_type = METHOD_TYPE_GFN1
      case ('gfn2', 'gfn2-xtb')
         method_type = METHOD_TYPE_GFN2

         ! SCF methods
      case ('hf', 'rhf', 'uhf', 'hartree-fock')
         method_type = METHOD_TYPE_HF
      case ('dft', 'ks', 'kohn-sham')
         method_type = METHOD_TYPE_DFT

         ! Multi-reference
      case ('mcscf', 'casscf', 'casci')
         method_type = METHOD_TYPE_MCSCF

         ! Perturbation theory
      case ('mp2', 'ri-mp2', 'df-mp2', 'scs-mp2', 'sos-mp2')
         method_type = METHOD_TYPE_MP2
      case ('mp2-f12', 'ri-mp2-f12', 'df-mp2-f12')
         method_type = METHOD_TYPE_MP2_F12

         ! Coupled cluster
      case ('ccsd', 'ri-ccsd', 'df-ccsd')
         method_type = METHOD_TYPE_CCSD
      case ('ccsd(t)', 'ri-ccsd(t)', 'df-ccsd(t)')
         method_type = METHOD_TYPE_CCSD_T
      case ('ccsd-f12', 'ri-ccsd-f12')
         method_type = METHOD_TYPE_CCSD_F12
      case ('ccsd(t)-f12', 'ri-ccsd(t)-f12')
         method_type = METHOD_TYPE_CCSD_T_F12

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
         ! Semi-empirical
      case (METHOD_TYPE_GFN1)
         method_str = "gfn1"
      case (METHOD_TYPE_GFN2)
         method_str = "gfn2"

         ! SCF methods
      case (METHOD_TYPE_HF)
         method_str = "hf"
      case (METHOD_TYPE_DFT)
         method_str = "dft"

         ! Multi-reference
      case (METHOD_TYPE_MCSCF)
         method_str = "mcscf"

         ! Perturbation theory
      case (METHOD_TYPE_MP2)
         method_str = "mp2"
      case (METHOD_TYPE_MP2_F12)
         method_str = "mp2-f12"

         ! Coupled cluster
      case (METHOD_TYPE_CCSD)
         method_str = "ccsd"
      case (METHOD_TYPE_CCSD_T)
         method_str = "ccsd(t)"
      case (METHOD_TYPE_CCSD_F12)
         method_str = "ccsd-f12"
      case (METHOD_TYPE_CCSD_T_F12)
         method_str = "ccsd(t)-f12"

      case default
         method_str = "unknown"
      end select

   end function method_type_to_string

end module mqc_method_types
