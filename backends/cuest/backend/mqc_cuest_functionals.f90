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
   use pic_ascii, only: to_lower
   implicit none
   private

   public :: functional_name_to_id   !! Functional name -> cuEST identifier
   public :: supported_functionals   !! Human-readable list, for error messages

contains

   subroutine vv10_unsupported(name, error)
      !! Refuse a `-V` functional on this backend
      !!
      !! **cuEST can do this and mqc does not ask it to.** The library exposes
      !! `CUEST_XCINTPLAN_IS_VV10`, the `VV10_B`/`VV10_C`/`VV10_SCALE`
      !! attributes and a separate `cuestNonlocalXCPotential*Compute` entry
      !! point; this backend sets the exchange scales and nothing else, and
      !! never calls that routine. The semilocal half alone converges and
      !! prints a number 43 mHa wrong on water -- so before this refusal, a GPU
      !! run of wb97x-v produced a confident wrong answer where the CPU path
      !! produced an error.
      !!
      !! Refused by name rather than by querying IS_VV10, which would be the
      !! better test and would cover functionals cuEST adds later, because it
      !! has to hold on a machine with no GPU to check it on. Wiring the
      !! nonlocal compute is the real fix and needs a card.
      character(len=*), intent(in) :: name
      type(error_t), intent(out) :: error

      call error%set(ERROR_VALIDATION, "'"//name//"' carries a non-local correlation "// &
                     "term (VV10). cuEST implements it, but this backend does not yet "// &
                     "enable it, and the semilocal part alone is a different functional "// &
                     "-- 43 mHa on water. Refused rather than evaluated without it; the "// &
                     "CPU path runs this functional.")
   end subroutine vv10_unsupported

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

      key = to_lower(adjustl(name))
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
         call vv10_unsupported("b97m-v", error)
         return
      case ("lc-wpbe", "lcwpbe")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_LCWPBE
      case ("wb97x")
         functional_id = CUEST_XCINTPLAN_PARAMETERS_FUNCTIONAL_WB97X
      case ("wb97x-v", "wb97xv")
         call vv10_unsupported("wb97x-v", error)
         return
      case ("wb97m-v", "wb97mv")
         call vv10_unsupported("wb97m-v", error)
         return
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
