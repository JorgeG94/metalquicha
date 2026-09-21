!! What a functional is called to the dispersion library, and which corrections exist
module mqc_dispersion_names
   !! Two spellings meet here and neither side may guess at the other.
   !!
   !! A deck names a functional in this program's spelling -- the one
   !! `xc_spec_from_name` parses -- and the dispersion library keys its damping
   !! parameters by
   !! its own. `camb3lyp` against `cam-b3lyp` is the harmless kind of
   !! difference; the dangerous kind is a functional this program supports that
   !! has no published D3 parametrisation at all, because handing the library a
   !! near neighbour's name returns a number rather than an error. Every such
   !! case is refused here by name, with the reason, and the list of what is
   !! known is part of the message.
   !!
   !! There is one table per correction and they are not shared. s-dftd3 and
   !! dftd4 are different libraries with different published parametrisations:
   !! a functional with D3 damping parameters need not have D4 ones, and
   !! reusing one table for the other would be the near-neighbour mistake in a
   !! new place.
   !!
   !! Separate from the wrappers that call the libraries so that it is compiled
   !! -- and tested -- whether or not this build has either of them. The mapping
   !! is a fact about two vocabularies, not about a link line.
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: d3_functional_alias
   public :: d4_functional_alias
   public :: dispersion_kind_is_known
   public :: dispersion_kind_is_d3
   public :: dispersion_kind_is_d4
   public :: DISPERSION_KINDS

   character(len=*), parameter :: DISPERSION_KINDS = "d3bj, d4"
      !! The values `keywords.dft.dispersion` accepts, besides `false`.
      !!
      !! Two, and each is served by a different library: "d3bj" by s-dftd3 and
      !! "d4" by dftd4, behind MQC_ENABLE_DFTD3 and MQC_ENABLE_DFTD4
      !! respectively. Zero damping is still reachable through s-dftd3's
      !! `dftd3_load_zero_damping` and belongs here when it is wired, not
      !! before: a keyword that accepts a word it then ignores is worse than one
      !! that refuses it.

   integer, parameter :: MAX_NAME = 32

   character(len=*), parameter :: KNOWN_LIST = &
                                  "b2gp-plyp, b2plyp, b3lyp, blyp, cam-b3lyp, m06-l, mpw2plyp, "// &
                                  "pbe, pbe0, r2scan, r2scan0, r2scan50, r2scanh, scan, tpss, wb97x"

   character(len=*), parameter :: D4_KNOWN_LIST = &
                                  "b2gp-plyp, b2plyp, b3lyp, blyp, cam-b3lyp, m06-l, mpw2plyp, "// &
                                  "pbe, pbe0, r2scan, r2scan0, r2scan50, r2scanh, scan, tpss, wb97x"
      !! Written out separately from `KNOWN_LIST` even where the two coincide.
      !!
      !! They coincide today over the functionals `xc_spec_from_name` accepts,
      !! which was checked against dftd4 v4.2.0's own tables rather than assumed
      !! -- every name below is exercised against the library in
      !! test/test_mqc_dispersion.f90. It is not a coincidence to rely on: dftd4
      !! carries D4 fits for wB97M and B97M that s-dftd3 has no D3 counterpart
      !! for, and s-dftd3 carries names dftd4 does not, so the day
      !! `xc_spec_from_name` learns one of them the two lists part company. One
      !! list shared between them would make that day silent.

contains

   pure function dispersion_kind_is_known(kind) result(known)
      !! Whether `keywords.dft.dispersion` named a correction this program has
      !!
      !! Spelling only. Whether the *build* linked the library that serves it is
      !! a different question, asked by `dispersion_kind_available` in
      !! `mqc_dispersion_apply`.
      character(len=*), intent(in) :: kind
      logical :: known

      known = dispersion_kind_is_d3(kind) .or. dispersion_kind_is_d4(kind)
   end function dispersion_kind_is_known

   pure function dispersion_kind_is_d3(kind) result(is_d3)
      !! Whether this correction is one s-dftd3 serves
      character(len=*), intent(in) :: kind
      logical :: is_d3

      is_d3 = (trim(adjustl(kind)) == "d3bj")
   end function dispersion_kind_is_d3

   pure function dispersion_kind_is_d4(kind) result(is_d4)
      !! Whether this correction is one dftd4 serves
      character(len=*), intent(in) :: kind
      logical :: is_d4

      is_d4 = (trim(adjustl(kind)) == "d4")
   end function dispersion_kind_is_d4

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

      alias = ""
      given = adjustl(functional)
      lower = folded(given)

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

   pure subroutine d4_functional_alias(functional, alias, error)
      !! dftd4's spelling of `functional`, or a refusal saying why there is none
      !!
      !! The same three kinds of refusal as `d3_functional_alias`, and a
      !! deliberately separate table: see the note on `D4_KNOWN_LIST`. A
      !! functional missing here because dftd4 never fitted it must reach the
      !! refusal, not the library, which answers a near neighbour's name with a
      !! number rather than an error.
      character(len=*), intent(in) :: functional
      character(len=MAX_NAME), intent(out) :: alias
      type(error_t), intent(out) :: error

      character(len=MAX_NAME) :: lower, given

      alias = ""
      given = adjustl(functional)
      lower = folded(given)

      select case (trim(lower))
         ! Spellings dftd4 shares with us, listed anyway rather than passed
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
         ! wB97X-D4, dftd4's own entry under this name. Not wB97X-V below, and
         ! not the same fit as s-dftd3's `wb97x` either -- which is exactly why
         ! the two tables are separate even where the spellings agree.
         alias = "wb97x"
      case ("b2plyp")
         alias = "b2plyp"
      case ("mpw2plyp", "mpw2-plyp")
         alias = "mpw2plyp"
         ! Where the two vocabularies differ. dftd4 accepts the hyphenated
         ! spellings too, but the unhyphenated ones are written here so that
         ! this table says what is sent rather than relying on the library's
         ! own aliasing.
      case ("cam-b3lyp", "camb3lyp")
         alias = "camb3lyp"
      case ("m06-l", "m06l")
         alias = "m06l"
      case ("b2gp-plyp", "b2gpplyp")
         alias = "b2gpplyp"
      case ("wb97x-v", "wb97xv", "wb97m-v", "wb97mv", "b97m-v", "b97mv")
         ! dftd4 does carry wB97M-D4 and B97M-D4, which s-dftd3 has no
         ! counterpart for -- but those are the functionals *without* VV10.
         ! `model.functional` names the -V forms, whose non-local correlation
         ! already accounts for dispersion, so the refusal is the same one.
         call error%set(ERROR_VALIDATION, "functional '"//trim(given)// &
                        "' already contains its own non-local correlation (VV10). Adding "// &
                        "an empirical D4 correction on top of it counts dispersion twice. "// &
                        "Drop keywords.dft.dispersion, or ask for a functional without "// &
                        "the -V.")
      case ("scan0", "r2scan01", "svwn", "svwn5", "lda", "lsda")
         call error%set(ERROR_VALIDATION, "functional '"//trim(given)// &
                        "' has no published D4 damping parameters, and a neighbouring "// &
                        "functional's would be a number with nothing behind it. "// &
                        "D4 is available for: "//D4_KNOWN_LIST//".")
      case default
         call error%set(ERROR_VALIDATION, "no D4 damping parameters are known here for "// &
                        "functional '"//trim(given)//"'. D4 is "// &
                        "available for: "//D4_KNOWN_LIST//".")
      end select
   end subroutine d4_functional_alias

   pure function folded(given) result(lower)
      !! `given` with A-Z folded to a-z, and nothing else touched
      !!
      !! A deck may shout, and both libraries key their tables in lower case.
      !! ASCII by hand rather than by intrinsic because Fortran has none.
      character(len=*), intent(in) :: given
      character(len=MAX_NAME) :: lower

      integer :: i, code

      lower = ""
      do i = 1, min(len_trim(given), MAX_NAME)
         code = iachar(given(i:i))
         if (code >= iachar("A") .and. code <= iachar("Z")) then
            lower(i:i) = achar(code + 32)
         else
            lower(i:i) = given(i:i)
         end if
      end do
   end function folded

end module mqc_dispersion_names
