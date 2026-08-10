!! Mapping of exchange-correlation functional names to cuEST identifiers
module mqc_cuest_functionals
   !! Translates the functional name from an input file into the identifier
   !! cuEST's XC integral plan expects.
   !!
   !! Only functionals cuEST implements natively are listed. The amount of
   !! exact exchange each one carries is deliberately *not* tabulated here:
   !! it is queried from the XC plan after construction, so the Coulomb side
   !! can never drift out of step with the functional definition.
   use mqc_error, only: error_t, ERROR_VALIDATION
   use cuest, only: CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_HF, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B3LYP1, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B3LYP5, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B97, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_BLYP, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M06L, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE0, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_R2SCAN, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_SVWN5, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B97MV, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_LCWPBE, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97X, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97XV, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97MV, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_LCWPBEH, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_CAMB3LYP, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_HSE06, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M06, &
                    CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M062X
   implicit none
   private

   public :: functional_name_to_id   !! Functional name -> cuEST identifier
   public :: supported_functionals   !! Human-readable list, for error messages

contains

   pure function lowercase(string) result(lower)
      !! Lowercase a string
      character(len=*), intent(in) :: string
      character(len=len(string)) :: lower
      integer :: i, code

      do i = 1, len(string)
         code = iachar(string(i:i))
         if (code >= iachar("A") .and. code <= iachar("Z")) then
            lower(i:i) = achar(code + 32)
         else
            lower(i:i) = string(i:i)
         end if
      end do
   end function lowercase

   subroutine functional_name_to_id(name, functional_id, error)
      !! Translate a functional name into the cuEST XC plan identifier
      !!
      !! Aliases follow common usage: bare "b3lyp" means the VWN5 variant
      !! that Gaussian popularized is *not* assumed -- cuEST exposes both, and
      !! "b3lyp" maps to B3LYP5 while "b3lyp1" selects the VWN1-RPA form.
      character(len=*), intent(in) :: name
      integer, intent(out) :: functional_id
      type(error_t), intent(out) :: error

      character(len=len(name)) :: key

      key = lowercase(adjustl(name))
      functional_id = -1

      select case (trim(key))
      case ("hf", "hartree-fock")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_HF
      case ("b3lyp", "b3lyp5")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B3LYP5
      case ("b3lyp1")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B3LYP1
      case ("b97")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B97
      case ("blyp")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_BLYP
      case ("m06-l", "m06l")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M06L
      case ("pbe")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE
      case ("pbe0")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_PBE0
      case ("r2scan")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_R2SCAN
      case ("svwn5", "lda", "svwn")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_SVWN5
      case ("b97m-v", "b97mv")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_B97MV
      case ("lc-wpbe", "lcwpbe")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_LCWPBE
      case ("wb97x")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97X
      case ("wb97x-v", "wb97xv")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97XV
      case ("wb97m-v", "wb97mv")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97MV
      case ("lc-wpbeh", "lcwpbeh")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_LCWPBEH
      case ("cam-b3lyp", "camb3lyp")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_CAMB3LYP
      case ("hse06")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_HSE06
      case ("m06")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M06
      case ("m06-2x", "m062x")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_M062X
      case default
         call error%set(ERROR_VALIDATION, "Unsupported functional '"//trim(name)// &
                        "'. cuEST provides: "//supported_functionals())
      end select
   end subroutine functional_name_to_id

   pure function supported_functionals() result(list)
      !! Comma-separated list of accepted functional names
      character(len=:), allocatable :: list

      list = "svwn5, blyp, pbe, m06-l, r2scan, b97, b3lyp, b3lyp1, pbe0, m06, "// &
             "m06-2x, hse06, cam-b3lyp, lc-wpbe, lc-wpbeh, wb97x, wb97x-v, "// &
             "b97m-v, wb97m-v, hf"
   end function supported_functionals

end module mqc_cuest_functionals
