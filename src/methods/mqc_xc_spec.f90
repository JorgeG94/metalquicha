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

   !> Components in the largest composition here. Double hybrids need four --
   !> exchange, correlation, and the two fractions that are not libxc's -- and
   !> nothing in view needs more.
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

   subroutine xc_spec_from_name(name, spec, error)
      !! A functional name to its composition
      !!
      !! Anything not named here is passed through as a single libxc component, on
      !! the assumption that libxc knows it. That is the right default: libxc
      !! carries several thousand functionals and this module should not shadow any
      !! of them. A name libxc also does not know fails at resolution, with libxc's
      !! own error, which is a better message than one invented here.
      character(len=*), intent(in) :: name
      type(xc_spec_t), intent(out) :: spec
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: lower
      integer :: i

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

         ! PBE is exchange plus correlation under one name libxc does not pair.
      case ("pbe")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("gga_x_pbe", 1.0_dp)
         spec%component(2) = xc_component_t("gga_c_pbe", 1.0_dp)

      case ("svwn", "lda", "lsda")
         spec%from_libxc = .false.
         spec%n_components = 2
         spec%component(1) = xc_component_t("lda_x", 1.0_dp)
         spec%component(2) = xc_component_t("lda_c_vwn", 1.0_dp)

         ! ---- everything else is libxc's ------------------------------------
      case default
         spec%from_libxc = .true.
         spec%n_components = 1
         spec%component(1) = xc_component_t(lower, 1.0_dp)
      end select
   end subroutine xc_spec_from_name

end module mqc_xc_spec
