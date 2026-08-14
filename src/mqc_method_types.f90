!! Method type constants for quantum chemistry methods
module mqc_method_types
   !! Defines integer constants for quantum chemistry methods to avoid string comparisons
   !! throughout the codebase. Provides conversion utilities between string
   !! representations and integer constants.
   use pic_types, only: int32, dp
   implicit none
   private

   ! Public constants - Semi-empirical
   public :: METHOD_TYPE_GFN1, METHOD_TYPE_GFN2
   ! Public constants - SCF methods
   public :: METHOD_TYPE_HF, METHOD_TYPE_DFT
   ! Public constants - Multi-reference
   public :: METHOD_TYPE_MCSCF
   ! Public constants - Classical, solving no wavefunction of its own
   public :: METHOD_TYPE_EFP2
   ! Public constants - Correlation methods
   public :: METHOD_TYPE_MP2, METHOD_TYPE_CCSD, METHOD_TYPE_CCSD_T
   public :: METHOD_TYPE_MP2_F12, METHOD_TYPE_CCSD_F12, METHOD_TYPE_CCSD_T_F12
   public :: METHOD_TYPE_UNKNOWN

   ! Public functions
   public :: method_type_from_string, method_type_to_string
   public :: parse_method_string  !! Input-file spelling -> method type
   public :: method_spin_scaling  !! Spin-component scaling a spelling implies
   public :: method_wants_density_fitting  !! Whether a spelling asks for RI

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

   !> Effective fragment potentials (60-69). Unlike every other entry here, this
   !> solves no wavefunction of its own: each fragment carries one already, computed
   !> when its potential was made, and what is evaluated is the interaction between
   !> them.
   integer(int32), parameter :: METHOD_TYPE_EFP2 = 60

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
         if (lower_str(i:i) >= "A" .and. lower_str(i:i) <= "Z") then
            lower_str(i:i) = achar(iachar(lower_str(i:i)) + 32)
         end if
      end do

      ! Match against known types
      select case (lower_str)
         ! Semi-empirical
      case ("gfn1", "gfn1-xtb")
         method_type = METHOD_TYPE_GFN1
      case ("gfn2", "gfn2-xtb")
         method_type = METHOD_TYPE_GFN2

         ! SCF methods
      case ("hf", "rhf", "uhf", "hartree-fock")
         method_type = METHOD_TYPE_HF
      case ("dft", "ks", "kohn-sham")
         method_type = METHOD_TYPE_DFT

         ! Multi-reference
      case ("mcscf", "casscf", "casci")
         method_type = METHOD_TYPE_MCSCF

         ! Effective fragment potentials. "efp" alone means EFP2: there is no EFP1
         ! here, and spelling it out is what a deck written for another program does.
      case ("efp2", "efp")
         method_type = METHOD_TYPE_EFP2

         ! Perturbation theory
      case ("mp2", "ri-mp2", "df-mp2", "scs-mp2", "sos-mp2")
         method_type = METHOD_TYPE_MP2
      case ("mp2-f12", "ri-mp2-f12", "df-mp2-f12")
         method_type = METHOD_TYPE_MP2_F12

         ! Coupled cluster
      case ("ccsd", "ri-ccsd", "df-ccsd")
         method_type = METHOD_TYPE_CCSD
      case ("ccsd(t)", "ri-ccsd(t)", "df-ccsd(t)")
         method_type = METHOD_TYPE_CCSD_T
      case ("ccsd-f12", "ri-ccsd-f12")
         method_type = METHOD_TYPE_CCSD_F12
      case ("ccsd(t)-f12", "ri-ccsd(t)-f12")
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

         ! Effective fragment potentials
      case (METHOD_TYPE_EFP2)
         method_str = "efp2"

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

   subroutine method_spin_scaling(method_str, use_scs, ss_scale, os_scale)
      !! The spin-component scaling a method name implies, if any
      !!
      !! "scs-mp2" and "sos-mp2" both parse to METHOD_TYPE_MP2, so the type
      !! alone cannot tell them from plain MP2 and a deck asking for either
      !! would run unscaled and report it as scaled. The spelling is the only
      !! place the intent survives, so it is read here, beside the parser that
      !! already has the lowercased string.
      !!
      !! SCS-MP2 is Grimme: 1.2 opposite-spin, one third same-spin. SOS-MP2 is
      !! Head-Gordon: 1.3 opposite-spin with the same-spin term dropped, which
      !! is what makes it cheaper rather than merely different.
      character(len=*), intent(in) :: method_str
      logical, intent(out) :: use_scs
      real(dp), intent(out) :: ss_scale, os_scale

      character(len=:), allocatable :: lower_str
      integer :: i

      lower_str = trim(adjustl(method_str))
      do i = 1, len(lower_str)
         if (lower_str(i:i) >= "A" .and. lower_str(i:i) <= "Z") then
            lower_str(i:i) = achar(iachar(lower_str(i:i)) + 32)
         end if
      end do

      use_scs = .false.
      ss_scale = 1.0_dp
      os_scale = 1.0_dp
      select case (lower_str)
      case ("scs-mp2")
         use_scs = .true.
         os_scale = 1.2_dp
         ss_scale = 1.0_dp/3.0_dp
      case ("sos-mp2")
         use_scs = .true.
         os_scale = 1.3_dp
         ss_scale = 0.0_dp
      case default
         ! Every other spelling, including plain "mp2", is unscaled. The values
         ! set above the select already say so; this is here to be explicit
         ! that no scaling is the deliberate answer rather than a gap.
      end select
   end subroutine method_spin_scaling

   function method_wants_density_fitting(method_str) result(wants_df)
      !! Whether a method name asks for the fitted route
      !!
      !! "ri-mp2" and "df-mp2" parse to METHOD_TYPE_MP2 exactly as "mp2" does,
      !! so without this a deck asking for either would run the conventional
      !! four-index transform and report it as RI -- right answer, wrong method,
      !! and no way to tell from the output. Same reason `method_spin_scaling`
      !! exists: the spelling carries intent the type discards.
      character(len=*), intent(in) :: method_str
      logical :: wants_df

      character(len=:), allocatable :: lower_str
      integer :: i

      lower_str = trim(adjustl(method_str))
      do i = 1, len(lower_str)
         if (lower_str(i:i) >= "A" .and. lower_str(i:i) <= "Z") then
            lower_str(i:i) = achar(iachar(lower_str(i:i)) + 32)
         end if
      end do

      wants_df = (index(lower_str, "ri-") == 1) .or. (index(lower_str, "df-") == 1)
   end function method_wants_density_fitting

   function parse_method_string(method_str) result(method_type)
      !! Parse method string from input file (e.g., "XTB-GFN1" -> gfn1)
      character(len=*), intent(in) :: method_str
      integer(int32) :: method_type

      character(len=:), allocatable :: lower_str, method_part
      integer :: dash_pos, i

      ! Convert to lowercase
      allocate (character(len=len_trim(method_str)) :: lower_str)
      lower_str = trim(adjustl(method_str))
      do i = 1, len(lower_str)
         if (lower_str(i:i) >= "A" .and. lower_str(i:i) <= "Z") then
            lower_str(i:i) = achar(iachar(lower_str(i:i)) + 32)
         end if
      end do

      ! Handle "XTB-GFN1" format -> extract "gfn1"
      if (index(lower_str, "xtb") > 0) then
         dash_pos = index(lower_str, "-")
         if (dash_pos > 0) then
            method_part = lower_str(dash_pos + 1:)
         else
            method_part = lower_str
         end if
      else
         method_part = lower_str
      end if

      method_type = method_type_from_string(method_part)

   end function parse_method_string

end module mqc_method_types
