!! What a functional name means, as a composition
module mqc_xc_spec
   !! Turns a functional name into the pieces a calculation has to assemble.
   !!
   !! For almost every functional this is a formality: libxc knows the name, owns
   !! the definition, and reports its own exact-exchange fraction through
   !! `xc_hyb_exx_coef`. Nothing here should duplicate any of that, and for those
   !! functionals a spec is one component with weight one.
   !!
   !! **Double hybrids are the exception, and the reason this module exists.**
   !! libxc 7.1.2 does not carry them -- not "exposes no coefficient for them", but
   !! carries no such functional at all: B2PLYP, DSD-PBEP86 and the XYG family are
   !! absent, and the only double-hybrid-looking entries are BHANDH and BHANDHLYP,
   !! which are plain fifty-percent hybrids. So a double hybrid cannot be asked
   !! for; it has to be *built*, out of libxc's individual exchange and
   !! correlation functionals, an exact-exchange fraction and a second-order
   !! perturbative fraction.
   !!
   !! Those weights are therefore **our data**, which is the whole risk: a table of
   !! coefficients drifts, gets copied wrong, and then produces a plausible number.
   !! Three things keep that honest, and all three are deliberate:
   !!
   !!   * every coefficient appears exactly once, here, with the paper it came
   !!     from written beside it;
   !!   * component functionals are named rather than numbered, so this table
   !!     carries no copy of libxc's own id list and a rename fails loudly at
   !!     resolution rather than silently selecting a different functional;
   !!   * the defining relation of a double hybrid -- that the semilocal exchange
   !!     and correlation weights are the complements of the exact-exchange and
   !!     perturbative fractions -- is asserted in the tests, so a mistyped row
   !!     does not survive.
   !!
   !! No libxc dependency. This module is a specification, not an evaluation, so it
   !! builds and is testable whether or not `MQC_ENABLE_LIBXC` is on. Resolving a
   !! component name to a libxc id, and evaluating it, belong to the layer that has
   !! libxc.
   use pic_types, only: dp
   use mqc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: xc_spec_t, xc_component_t
   public :: xc_spec_from_name
   public :: MAX_XC_COMPONENTS

   !! Components in the largest composition here. Double hybrids need four --
   !! exchange, correlation, and the two fractions that are not libxc's -- and
   !! nothing in view needs more.
   integer, parameter :: MAX_XC_COMPONENTS = 6

   type :: xc_component_t
      !! One libxc functional and how much of it to take
      character(len=48) :: name = ""
         !! libxc's own name, e.g. "gga_x_b88". Resolved to an id by the layer
         !! that has libxc, so this table never holds a number that could drift
         !! out of step with the library.
      real(dp) :: weight = 1.0_dp
   end type xc_component_t

   type :: xc_spec_t
      !! A functional, as the sum of things that can actually be evaluated
      !!
      !!     E_xc = sum_i weight_i * E[component_i]
      !!            + exx_fraction * E_x^exact
      !!            + pt2_fraction * E_c^MP2
      !!
      !! `exx_fraction` is left at zero for anything libxc will report itself:
      !! for a plain hybrid the evaluation layer should ask `xc_hyb_exx_coef`
      !! rather than trust a number here, and a nonzero value means the name is
      !! one libxc does not know.
      character(len=32) :: name = ""
      integer :: n_components = 0
      type(xc_component_t) :: component(MAX_XC_COMPONENTS)
      real(dp) :: exx_fraction = 0.0_dp
      real(dp) :: pt2_fraction = 0.0_dp
      logical :: from_libxc = .true.
         !! True when libxc owns the whole definition and this spec is only
         !! passing the name through. False for a composition assembled here,
         !! which is the case that needs a citation and a test.
   contains
      procedure :: is_double_hybrid => spec_is_double_hybrid
      procedure :: needs_exact_exchange => spec_needs_exact_exchange
   end type xc_spec_t

contains

   pure function spec_is_double_hybrid(this) result(is_dh)
      !! Whether this functional needs an MP2 correlation term
      class(xc_spec_t), intent(in) :: this
      logical :: is_dh
      is_dh = this%pt2_fraction /= 0.0_dp
   end function spec_is_double_hybrid

   pure function spec_needs_exact_exchange(this) result(needs)
      !! Whether a Fock exchange build is required
      !!
      !! True for a composition that states its own fraction. For a libxc-owned
      !! hybrid this is false and the evaluation layer must ask libxc, because the
      !! fraction is libxc's to report and copying it here would be a second
      !! source of truth.
      class(xc_spec_t), intent(in) :: this
      logical :: needs
      needs = this%exx_fraction /= 0.0_dp
   end function spec_needs_exact_exchange

   subroutine xc_spec_from_name(name, spec, error, allow_half)
      !! A functional name to its composition
      !!
      !! Anything not named here is passed through as a single libxc component, on
      !! the assumption that libxc knows it. That is the right default: libxc
      !! carries several thousand functionals and this module should not shadow any
      !! of them. A name libxc also does not know fails at resolution, with libxc's
      !! own error, which is a better message than one invented here.
      !!
      !! The one thing checked before that pass-through is that the name is a
      !! *whole* functional and not one half of one -- see `check_name_is_whole`,
      !! which is the difference between a deck being refused and a deck quietly
      !! returning the wrong energy.
      character(len=*), intent(in) :: name
      type(xc_spec_t), intent(out) :: spec
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: allow_half
         !! Accept a libxc name that is only exchange or only correlation.
         !! Default false, and a deck must never set it: the point of the check
         !! is that half a functional converges silently. It exists for the
         !! derivative tests, which want one rung's kernel in isolation --
         !! `lda_x` before `svwn`, `gga_x_pbe` before `pbe` -- so that a broken
         !! exchange second derivative is not masked by a working correlation
         !! one. Naming it at the call site is the whole value: a half asked for
         !! on purpose reads differently from a half asked for by mistake.

      character(len=:), allocatable :: lower
      logical :: whole_only
      integer :: i

      whole_only = .true.
      if (present(allow_half)) whole_only = .not. allow_half

      lower = trim(adjustl(name))
      do i = 1, len(lower)
         if (lower(i:i) >= "A" .and. lower(i:i) <= "Z") then
            lower(i:i) = achar(iachar(lower(i:i)) + 32)
         end if
      end do

      if (len(lower) == 0) then
         call error%set(ERROR_VALIDATION, "no exchange-correlation functional was named")
         return
      end if

      spec%name = lower

      select case (lower)

         ! ---- double hybrids, which libxc does not carry ---------------------
         !
         ! B2PLYP: Grimme, J. Chem. Phys. 124, 034108 (2006). Half-and-half in
         ! spirit: a fraction a_x of exact exchange with the complement from B88,
         ! and a fraction a_c of MP2 correlation with the complement from LYP.
      case ("b2plyp")
         spec%from_libxc = .false.
         spec%exx_fraction = 0.53_dp
         spec%pt2_fraction = 0.27_dp
         spec%n_components = 2
         spec%component(1) = xc_component_t("gga_x_b88", 1.0_dp - spec%exx_fraction)
         spec%component(2) = xc_component_t("gga_c_lyp", 1.0_dp - spec%pt2_fraction)

         ! B2GP-PLYP: Karton et al., J. Phys. Chem. A 112, 12868 (2008).
      case ("b2gp-plyp", "b2gpplyp")
         spec%from_libxc = .false.
         spec%exx_fraction = 0.65_dp
         spec%pt2_fraction = 0.36_dp
         spec%n_components = 2
         spec%component(1) = xc_component_t("gga_x_b88", 1.0_dp - spec%exx_fraction)
         spec%component(2) = xc_component_t("gga_c_lyp", 1.0_dp - spec%pt2_fraction)

         ! mPW2-PLYP: Schwabe and Grimme, PCCP 8, 4398 (2006).
      case ("mpw2plyp", "mpw2-plyp")
         spec%from_libxc = .false.
         spec%exx_fraction = 0.55_dp
         spec%pt2_fraction = 0.25_dp
         spec%n_components = 2
         spec%component(1) = xc_component_t("gga_x_mpw91", 1.0_dp - spec%exx_fraction)
         spec%component(2) = xc_component_t("gga_c_lyp", 1.0_dp - spec%pt2_fraction)

         ! ---- compositions libxc has the parts for but not the name ----------
         !
         ! SVWN is the textbook local functional -- Slater exchange with VWN
         ! correlation -- and libxc carries both pieces without carrying that
         ! name for the pair. Listed here rather than made a special case in the
         ! evaluator, because "a name libxc lacks, built from parts it has" is
         ! exactly what this table is for; double hybrids are only the hardest
         ! instance of it.
         ! Hybrid aliases. libxc spells these hyb_gga_xc_b3lyp and hyb_gga_xc_pbeh,
         ! and carries the exact-exchange fraction itself -- so these rows are pure
         ! renames and must state no fraction of their own, or there would be two
         ! sources of truth for it.
      case ("b3lyp")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_gga_xc_b3lyp", 1.0_dp)

         ! Range-separated hybrids, by their libxc names.
      case ("wb97x")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_gga_xc_wb97x", 1.0_dp)

         ! The -V family, which libxc carries and this table used to omit. Named
         ! here so an unsupported one is refused for the reason it is actually
         ! unsupported -- the non-local correlation check in `mqc_libcint_xc`, or
         ! the meta-GGA one -- rather than with "libxc does not know the
         ! functional", which was never true and read like a typo in the deck.
      case ("wb97x-v", "wb97xv")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_gga_xc_wb97x_v", 1.0_dp)

      case ("wb97m-v", "wb97mv")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_mgga_xc_wb97m_v", 1.0_dp)

      case ("b97m-v", "b97mv")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("mgga_xc_b97m_v", 1.0_dp)

      case ("cam-b3lyp", "camb3lyp")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_gga_xc_cam_b3lyp", 1.0_dp)

      case ("pbe0", "pbeh")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_gga_xc_pbeh", 1.0_dp)

         ! Meta-GGAs, same pattern: libxc carries the halves, not the pair.
      case ("tpss")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("mgga_x_tpss", 1.0_dp)
         spec%component(2) = xc_component_t("mgga_c_tpss", 1.0_dp)

      case ("m06-l", "m06l")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("mgga_x_m06_l", 1.0_dp)
         spec%component(2) = xc_component_t("mgga_c_m06_l", 1.0_dp)

         ! ---- the SCAN family ------------------------------------------------
         !
         ! SCAN: Sun, Ruzsinszky and Perdew, Phys. Rev. Lett. 115, 036402 (2015).
         ! r2SCAN: Furness, Kaplan, Ning, Perdew and Sun, J. Phys. Chem. Lett. 11,
         ! 8208 (2020) -- SCAN re-regularised to remove the numerical sensitivity
         ! that made the original need an unusually fine grid. It is the member of
         ! the family to reach for, and the reason the rest are spelled out beside
         ! it rather than left to whoever knows libxc's naming.
         !
         ! The semilocal members follow TPSS and M06-L above: libxc carries the two
         ! halves and no name for the pair. Until these rows existed `"functional":
         ! "r2scan"` failed as a name libxc does not know -- while `mgga_x_r2scan`
         ! *succeeded*, converged, at an energy missing every bit of the
         ! correlation. That second failure is the one worth preventing, and the
         ! completeness check in the default case below is what now prevents it.
      case ("r2scan")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("mgga_x_r2scan", 1.0_dp)
         spec%component(2) = xc_component_t("mgga_c_r2scan", 1.0_dp)

         ! The same construction at the larger regularisation parameter. libxc
         ! labels both halves "with larger value for eta", which is what says they
         ! are a pair -- there is no published name pairing them for us.
      case ("r2scan01")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("mgga_x_r2scan01", 1.0_dp)
         spec%component(2) = xc_component_t("mgga_c_r2scan01", 1.0_dp)

      case ("scan")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("mgga_x_scan", 1.0_dp)
         spec%component(2) = xc_component_t("mgga_c_scan", 1.0_dp)

         ! The hybrids split two ways, and the split is libxc's rather than ours.
         ! r2SCAN0, r2SCANh and r2SCAN50 it carries whole, one `hyb_mgga_xc_` name
         ! each, so those three rows are renames and must state no fraction of
         ! their own. SCAN0 it does not: `hyb_mgga_x_scan0` is the hybridised
         ! exchange by itself, so SCAN0 is that plus SCAN correlation -- and the
         ! quarter of exact exchange still belongs to libxc, which reports it from
         ! the component, not from a number written here.
      case ("scan0")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("hyb_mgga_x_scan0", 1.0_dp)
         spec%component(2) = xc_component_t("mgga_c_scan", 1.0_dp)

      case ("r2scan0")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_mgga_xc_r2scan0", 1.0_dp)

      case ("r2scanh")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_mgga_xc_r2scanh", 1.0_dp)

      case ("r2scan50")
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t("hyb_mgga_xc_r2scan50", 1.0_dp)

         ! BLYP, which libxc also leaves as two halves. Worth a row of its own
         ! rather than leaving callers to spell both: it is B2PLYP's semilocal
         ! part at full weight, so it is the functional a double hybrid's kernel
         ! and gradient are checked against.
      case ("blyp")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("gga_x_b88", 1.0_dp)
         spec%component(2) = xc_component_t("gga_c_lyp", 1.0_dp)

         ! PBE is exchange plus correlation under one name libxc does not pair.
      case ("pbe")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("gga_x_pbe", 1.0_dp)
         spec%component(2) = xc_component_t("gga_c_pbe", 1.0_dp)

      case ("svwn", "svwn5", "lda", "lsda")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("lda_x", 1.0_dp)
         spec%component(2) = xc_component_t("lda_c_vwn", 1.0_dp)

         ! ---- everything else is libxc's ------------------------------------
      case default
         if (whole_only) call check_name_is_whole(lower, error)
         if (error%has_error()) return
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t(lower, 1.0_dp)
      end select
   end subroutine xc_spec_from_name

   pure function libxc_role(name) result(role)
      !! Which part of a functional a libxc name is: "x", "c", "xc", "k" or ""
      !!
      !! libxc spells every functional `<family>_<role>_<label>`, and the role is
      !! the first underscore-delimited segment that is one of those four. All 702
      !! functionals in 7.1.2 follow it without exception, so this reads libxc's
      !! own convention rather than keeping a list that could fall out of step
      !! with it. A name following no such convention returns "" and is left for
      !! libxc to reject in its own words.
      character(len=*), intent(in) :: name
      character(len=2) :: role

      integer :: i, lo, n

      role = ""
      lo = 1
      n = len_trim(name)
      do i = 1, n + 1
         if (i <= n) then
            if (name(i:i) /= "_") cycle
         end if
         if (is_role(name(lo:i - 1))) then
            role = name(lo:i - 1)
            return
         end if
         lo = i + 1
      end do
   end function libxc_role

   pure function is_role(segment) result(yes)
      !! Whether one name segment is a role marker rather than part of the label
      character(len=*), intent(in) :: segment
      logical :: yes
      yes = segment == "x" .or. segment == "c" .or. segment == "xc" .or. segment == "k"
   end function is_role

   subroutine check_name_is_whole(name, error)
      !! Refuse a libxc name that is only half of a functional
      !!
      !! libxc splits most of its semilocal functionals into an exchange half and
      !! a correlation half, and either half alone is a perfectly valid libxc
      !! functional that initialises, evaluates, and converges. It just is not the
      !! functional anyone meant: `mgga_x_r2scan` alone lands 0.32 Hartree from
      !! r2SCAN and `mgga_c_r2scan` alone nine Hartree, both silently, on a water
      !! molecule. Nothing downstream can catch that, because there is nothing
      !! wrong with the calculation -- only with which functional it ran.
      !!
      !! This costs nothing in reach: a pair libxc does not carry under one name
      !! belongs in the table above, which is where every such functional in this
      !! module already is.
      character(len=*), intent(in) :: name
      type(error_t), intent(inout) :: error

      character(len=2) :: role

      role = libxc_role(name)
      select case (trim(role))
      case ("x", "c")
         call error%set(ERROR_VALIDATION, "'"//trim(name)//"' is only the "// &
                        trim(merge("exchange   ", "correlation", trim(role) == "x"))// &
                        " half of a functional. libxc carries the halves "// &
                        "separately and one of them alone converges to an energy "// &
                        "that is not the functional it is named after. Ask for a "// &
                        "name that pairs them, or for a libxc name carrying "// &
                        "'_xc_'.")
      case ("k")
         call error%set(ERROR_VALIDATION, "'"//trim(name)//"' is a kinetic energy "// &
                        "functional, not an exchange-correlation one.")
      case default
         ! "xc", a whole functional; or "", a name following no libxc convention
         ! at all. Both pass, the second because an unrecognisable name is libxc's
         ! to reject in its own words rather than this module's to guess at.
      end select
   end subroutine check_name_is_whole

end module mqc_xc_spec
