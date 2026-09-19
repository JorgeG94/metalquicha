!! What a functional is called to the dispersion library, and which corrections exist
module mqc_dispersion_names
   !! Two spellings meet here and neither side may guess at the other.
   !!
   !! A deck names a functional in this program's spelling -- the one
   !! `xc_spec_from_name` parses -- and s-dftd3 keys its damping parameters by
   !! its own. `camb3lyp` against `cam-b3lyp` is the harmless kind of
   !! difference; the dangerous kind is a functional this program supports that
   !! has no published D3 parametrisation at all, because handing the library a
   !! near neighbour's name returns a number rather than an error. Every such
   !! case is refused here by name, with the reason, and the list of what is
   !! known is part of the message.
   !!
   !! Separate from the wrapper that calls the library so that it is compiled --
   !! and tested -- whether or not this build has s-dftd3. The mapping is a fact
   !! about two vocabularies, not about a link line.
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: d3_functional_alias
   public :: dispersion_kind_is_known
   public :: DISPERSION_KINDS

   character(len=*), parameter :: DISPERSION_KINDS = "d3bj"
      !! The values `keywords.dft.dispersion` accepts, besides `false`.
      !!
      !! One, for now. Zero damping and D4 are both reachable through the same
      !! library -- `dftd3_load_zero_damping`, and dftd4's own C API -- and
      !! belong here when they are wired, not before: a keyword that accepts a
      !! word it then ignores is worse than one that refuses it.

   integer, parameter :: MAX_NAME = 32

   character(len=*), parameter :: KNOWN_LIST = &
                                  "b2gp-plyp, b2plyp, b3lyp, blyp, cam-b3lyp, m06-l, mpw2plyp, "// &
                                  "pbe, pbe0, r2scan, r2scan0, r2scan50, r2scanh, scan, tpss, wb97x"

contains

   pure function dispersion_kind_is_known(kind) result(known)
      !! Whether `keywords.dft.dispersion` named a correction this build can run
      character(len=*), intent(in) :: kind
      logical :: known

      known = (trim(adjustl(kind)) == "d3bj")
   end function dispersion_kind_is_known

   pure subroutine d3_functional_alias(functional, alias, error)
      !! s-dftd3's spelling of `functional`, or a refusal saying why there is none
      !!
      !! The refusals are the point of the routine. Three kinds are distinguished
      !! because they mean different things to whoever wrote the deck:
      !!
      !! * a functional whose -V variant already carries non-local correlation,
      !!   where adding D3 would count dispersion twice;
      !! * a functional with no published D3 parametrisation, where there is
      !!   simply no number to use;
      !! * a functional this program does not know at all, which is a typo.
      character(len=*), intent(in) :: functional
      character(len=MAX_NAME), intent(out) :: alias
      type(error_t), intent(out) :: error

      character(len=MAX_NAME) :: lower, given
      integer :: i, code

      alias = ""
      given = adjustl(functional)
      lower = ""
      do i = 1, len_trim(given)
         code = iachar(given(i:i))
         if (code >= iachar("A") .and. code <= iachar("Z")) then
            lower(i:i) = achar(code + 32)
         else
            lower(i:i) = given(i:i)
         end if
      end do

      select case (trim(lower))
         ! Spellings s-dftd3 shares with us, listed anyway rather than passed
         ! through: an unlisted name must reach the refusal below, not the
         ! library, so that "known to mqc" and "known here" cannot drift.
      case ("b3lyp")
         alias = "b3lyp"
      case ("blyp")
         alias = "blyp"
      case ("pbe")
         alias = "pbe"
      case ("pbe0", "pbeh")
         alias = "pbe0"
      case ("tpss")
         alias = "tpss"
      case ("scan")
         alias = "scan"
      case ("r2scan")
         alias = "r2scan"
      case ("r2scan0")
         alias = "r2scan0"
      case ("r2scanh")
         alias = "r2scanh"
      case ("r2scan50")
         alias = "r2scan50"
      case ("wb97x")
         ! wB97X-D3(BJ), Najibi and Goerigk's reparametrisation, which is what
         ! s-dftd3 carries under this name. Not the same as wB97X-V below.
         alias = "wb97x"
      case ("b2plyp")
         alias = "b2plyp"
      case ("mpw2plyp", "mpw2-plyp")
         alias = "mpw2plyp"
         ! Where the two vocabularies differ: hyphens s-dftd3 does not use.
      case ("cam-b3lyp", "camb3lyp")
         alias = "camb3lyp"
      case ("m06-l", "m06l")
         alias = "m06l"
      case ("b2gp-plyp", "b2gpplyp")
         alias = "b2gpplyp"
      case ("wb97x-v", "wb97xv", "wb97m-v", "wb97mv", "b97m-v", "b97mv")
         call error%set(ERROR_VALIDATION, "functional '"//trim(given)// &
                        "' already contains its own non-local correlation (VV10). Adding "// &
                        "an empirical D3 correction on top of it counts dispersion twice. "// &
                        "Drop keywords.dft.dispersion, or ask for a functional without "// &
                        "the -V.")
      case ("scan0", "r2scan01", "svwn", "svwn5", "lda", "lsda")
         call error%set(ERROR_VALIDATION, "functional '"//trim(given)// &
                        "' has no published D3 damping parameters, and a neighbouring "// &
                        "functional's would be a number with nothing behind it. "// &
                        "Dispersion is available for: "//KNOWN_LIST//".")
      case default
         call error%set(ERROR_VALIDATION, "no D3 damping parameters are known here for "// &
                        "functional '"//trim(given)//"'. Dispersion is "// &
                        "available for: "//KNOWN_LIST//".")
      end select
   end subroutine d3_functional_alias

end module mqc_dispersion_names
