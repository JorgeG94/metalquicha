!! What a fragmentation scheme may hand a method, and the one refusal site
module mqc_fragment_capabilities
   !! FMO, EE-MBE and EFMO solve their fragments inside the backend rather
   !! than through `qc_method_t`, so what a method may do there is a property
   !! of the pairing of scheme and method. `fragment_capabilities` states it
   !! and `fragment_refusal` is the one place a deck asking for more is turned
   !! away. See `mqc_docs/source/developer_fragment_solver.rst`.
   use mqc_calc_types, only: CALC_TYPE_ENERGY, calc_type_to_string
   use mqc_method_config, only: method_config_t
   use mqc_method_types, only: METHOD_TYPE_HF, METHOD_TYPE_DFT, METHOD_TYPE_MP2, &
                               method_type_to_string
   use mqc_error, only: error_t
   use mqc_xc_spec, only: xc_spec_t, xc_spec_from_name
   implicit none
   private

   integer, parameter, public :: FRAGMENT_SCHEME_FMO = 1
   integer, parameter, public :: FRAGMENT_SCHEME_EE_MBE = 2
   integer, parameter, public :: FRAGMENT_SCHEME_EFMO = 3

   public :: fragment_needs_t
   public :: fragment_capabilities_t
   public :: fragment_capabilities
   public :: fragment_refusal

   type :: fragment_needs_t
      !! What a deck asks of the fragment solver
      logical :: cut = .false.
         !! `keywords.fragmentation.bond_breaking = "afo"`: a Fock projector
         !! at a detached bond
      logical :: unrestricted = .false.
         !! `model.unrestricted`
      integer :: calc_type = CALC_TYPE_ENERGY
         !! The deck's driver, as a `CALC_TYPE_*` constant
      logical :: pieda = .false.
         !! `keywords.fragmentation.pieda`
      logical :: pieda_dispersion = .false.
         !! `keywords.fragmentation.pieda_dispersion` is not `"none"`
      logical :: dispersion = .false.
         !! `keywords.dft.dispersion` (`method_config%dft%use_dispersion`).
         !! Meaningless off a Kohn-Sham method, where it is always false.
   end type fragment_needs_t

   type :: fragment_capabilities_t
      !! What a method can do inside one fragmentation scheme
      !!
      !! Every method that runs at all accepts an embedding operator, so that
      !! is not a field.
      logical :: runs = .false.
         !! The method runs under this scheme. The other fields are false
         !! when this is.
      logical :: cut = .false.
      logical :: unrestricted = .false.
      logical :: gradient = .false.
      logical :: pieda = .false.
         !! The method's own PIEDA terms exist. Which schemes offer PIEDA at
         !! all is decided by `check_pieda_support` in `mqc_config_adapter`.
      logical :: pieda_dispersion = .false.
         !! PIEDA's empirical `Edi` can be added beside the method's terms.
         !! False for an MP2-family method, whose `Edi` is its correlation and
         !! would count dispersion twice.
      logical :: dispersion = .false.
         !! Empirical dispersion runs alongside this method under this
         !! scheme, added to every fragment and n-mer as GAMESS does. True for
         !! Kohn-Sham under FMO and EE-MBE; false for Hartree-Fock, where there
         !! is no functional to take damping parameters from, for the MP2
         !! family, and for EFMO.
   end type fragment_capabilities_t

contains

   function fragment_capabilities(scheme, config) result(cap)
      !! What `config%method_type` may do under `scheme`
      !!
      !! FMO and EE-MBE: Hartree-Fock or Kohn-Sham, with or without a detached
      !! bond, plus MP2, SCS/SOS-MP2 or RI-MP2 as correlation on the embedded
      !! Hartree-Fock reference, without a detached bond (the frozen orbitals
      !! at a cut would be correlated). A double hybrid does not run: its PT2
      !! part is not added. All of them carry PIEDA. Empirical dispersion
      !! (`keywords.dft.dispersion`) runs with Kohn-Sham only, per fragment and
      !! n-mer. PIEDA's own empirical `Edi` is offered with Hartree-Fock and
      !! Kohn-Sham, and not with the MP2 family.
      !! EFMO: Hartree-Fock, and MP2 or RI-MP2 without a detached bond (the
      !! frozen virtual at a cut would be correlated) and without spin-component
      !! scaling. Nothing runs unrestricted or with a gradient.
      integer, intent(in) :: scheme
      type(method_config_t), intent(in) :: config
      type(fragment_capabilities_t) :: cap

      select case (scheme)
      case (FRAGMENT_SCHEME_FMO, FRAGMENT_SCHEME_EE_MBE)
         if (config%method_type == METHOD_TYPE_HF) then
            cap%runs = .true.
            cap%cut = .true.
            cap%pieda = .true.
            cap%pieda_dispersion = .true.
         else if (config%method_type == METHOD_TYPE_DFT) then
            cap%runs = .not. double_hybrid(config%dft%functional)
            cap%cut = cap%runs
            cap%pieda = cap%runs
            cap%pieda_dispersion = cap%runs
            cap%dispersion = cap%runs
         else if (config%method_type == METHOD_TYPE_MP2) then
            cap%runs = .true.
            cap%pieda = .true.
         end if
      case (FRAGMENT_SCHEME_EFMO)
         select case (config%method_type)
         case (METHOD_TYPE_HF)
            cap%runs = .true.
            cap%cut = .true.
            cap%pieda = .true.
            cap%pieda_dispersion = .true.
         case (METHOD_TYPE_MP2)
            if (.not. config%corr%use_scs) cap%runs = .true.
         case default
            ! Nothing else runs.
         end select
      case default
         ! An unknown scheme runs nothing.
      end select
   end function fragment_capabilities

   function fragment_refusal(scheme, config, needs) result(message)
      !! Why this deck cannot run under `scheme`, or "" when it can
      !!
      !! Reports the first unmet need, checked in this order: the method at
      !! all, a non-Energy driver, `unrestricted`, a detached bond,
      !! dispersion, PIEDA, PIEDA's empirical dispersion, and dispersion
      !! together with PIEDA's.
      !!
      !! Whether the build has the library that serves the correction is not
      !! asked here: the deck reader refuses it by name, with the CMake
      !! option, before any scheme is chosen.
      integer, intent(in) :: scheme
      type(method_config_t), intent(in) :: config
      type(fragment_needs_t), intent(in) :: needs
      character(len=:), allocatable :: message

      type(fragment_capabilities_t) :: cap

      cap = fragment_capabilities(scheme, config)
      message = ""

      if (.not. cap%runs) then
         message = runs_refusal(scheme, config)
      else if (needs%calc_type /= CALC_TYPE_ENERGY .and. .not. cap%gradient) then
         message = driver_refusal(scheme, needs%calc_type)
      else if (needs%unrestricted .and. .not. cap%unrestricted) then
         message = unrestricted_refusal(scheme)
      else if (needs%cut .and. .not. cap%cut) then
         message = cut_refusal(scheme, config)
      else if (needs%dispersion .and. .not. cap%dispersion) then
         message = dispersion_refusal(scheme, config%method_type)
      else if (needs%pieda .and. .not. cap%pieda) then
         message = "keywords.fragmentation.pieda has no decomposition for "// &
                   "model.method '"//trim(method_type_to_string(config%method_type))// &
                   "' yet. Drop pieda, or set model.method to 'hf'."
      else if (needs%pieda_dispersion .and. .not. cap%pieda_dispersion) then
         message = "keywords.fragmentation.pieda_dispersion adds empirical "// &
                   "dispersion to a pair, and under model.method '"// &
                   trim(method_type_to_string(config%method_type))// &
                   "' the pair's Edi is its correlation energy, which already holds "// &
                   "the dispersion: both would count it twice. Set "// &
                   "keywords.fragmentation.pieda_dispersion to 'none'."
      else if (needs%pieda_dispersion .and. needs%dispersion) then
         message = "keywords.dft.dispersion adds the empirical dispersion to every "// &
                   "fragment and n-mer, so the dispersion interaction is already in "// &
                   "each pair's energy, and keywords.fragmentation.pieda_dispersion "// &
                   "would count it a second time. Set "// &
                   "keywords.fragmentation.pieda_dispersion to 'none'; PIEDA then "// &
                   "reports the dispersion as its Edi term."
      end if
   end function fragment_refusal

   function runs_refusal(scheme, config) result(message)
      !! Why `config%method_type` does not run under `scheme` at all
      integer, intent(in) :: scheme
      type(method_config_t), intent(in) :: config
      character(len=:), allocatable :: message

      if (scheme == FRAGMENT_SCHEME_EFMO) then
         if (config%method_type == METHOD_TYPE_MP2 .and. config%corr%use_scs) then
            message = "EFMO does not scale the spin components of its fragment MP2 "// &
                      "energies. SCS-MP2 fragments would be a different method from "// &
                      "the one the paper runs, so it is refused rather than quietly "// &
                      "given plain MP2."
         else
            message = "EFMO runs Hartree-Fock, MP2 or RI-MP2 fragments, and "// &
                      "model.method is '"// &
                      trim(method_type_to_string(config%method_type))// &
                      "'. Kohn-Sham fragments would need a MAKEFP that is not "// &
                      "restricted to a Hartree-Fock reference, and coupled-cluster "// &
                      "fragments are not implemented."
         end if
      else if (config%method_type == METHOD_TYPE_DFT) then
         message = "model.functional '"//trim(config%dft%functional)//"' is a double "// &
                   "hybrid, and FMO and EE-MBE do not add its perturbative "// &
                   "correlation. Choose a functional with no PT2 part."
      else
         message = "The fragment calculations of FMO and EE-MBE currently run "// &
                   "Hartree-Fock, Kohn-Sham, and MP2, SCS-MP2, SOS-MP2 or RI-MP2 as "// &
                   "correlation on top, and model.method is '"// &
                   trim(method_type_to_string(config%method_type))//"', which is "// &
                   "not yet wired into them. Set model.method to 'hf', a functional, "// &
                   "or 'mp2'/'ri-mp2'."
      end if
   end function runs_refusal

   function driver_refusal(scheme, calc_type) result(message)
      !! Why a non-Energy driver is refused under `scheme`
      integer, intent(in) :: scheme
      integer, intent(in) :: calc_type
      character(len=:), allocatable :: message

      if (scheme == FRAGMENT_SCHEME_EFMO) then
         message = "EFMO computes energies only, and driver is '"// &
                   trim(calc_type_to_string(calc_type))//"'. EFMO gradients are "// &
                   "not implemented, with or without detached bonds."
      else
         message = "FMO and EE-MBE compute energies only, and driver is '"// &
                   trim(calc_type_to_string(calc_type))//"'. Fragment gradients "// &
                   "and Hessians are not implemented; set driver to 'Energy'."
      end if
   end function driver_refusal

   function unrestricted_refusal(scheme) result(message)
      !! Why `model.unrestricted` is refused under `scheme`
      integer, intent(in) :: scheme
      character(len=:), allocatable :: message

      if (scheme == FRAGMENT_SCHEME_EFMO) then
         message = "EFMO is closed-shell for now: every fragment and every dimer "// &
                   "is solved with restricted Hartree-Fock, so "// &
                   "model.unrestricted cannot be honoured."
      else
         message = "FMO and EE-MBE are closed-shell for now: every fragment and "// &
                   "every n-mer is solved restricted, so "// &
                   "model.unrestricted cannot be honoured."
      end if
   end function unrestricted_refusal

   function double_hybrid(functional) result(is_dh)
      !! Whether `functional` names a double hybrid
      !!
      !! A name `xc_spec_from_name` rejects is not a double hybrid here; the
      !! backend reports it when it builds the functional.
      character(len=*), intent(in) :: functional
      logical :: is_dh

      type(xc_spec_t) :: spec
      type(error_t) :: error

      call xc_spec_from_name(functional, spec, error)
      is_dh = .false.
      if (.not. error%has_error()) is_dh = spec%is_double_hybrid()
   end function double_hybrid

   function dispersion_refusal(scheme, method_type) result(message)
      !! Why `keywords.dft.dispersion` is refused under `scheme`
      integer, intent(in) :: scheme
      integer, intent(in) :: method_type
      character(len=:), allocatable :: message

      if (scheme == FRAGMENT_SCHEME_EFMO) then
         message = "EFMO does not add empirical dispersion: its fragments are "// &
                   "Hartree-Fock or MP2, and dispersion is added to a Kohn-Sham "// &
                   "reference. Drop keywords.dft.dispersion."
      else if (method_type == METHOD_TYPE_MP2) then
         message = trim(scheme_name(scheme))//" adds empirical dispersion to "// &
                   "Kohn-Sham fragments and n-mers, and model.method '"// &
                   trim(method_type_to_string(method_type))//"' already holds the "// &
                   "dispersion in its correlation energy. Drop keywords.dft.dispersion."
      else
         message = trim(scheme_name(scheme))//" adds empirical dispersion to "// &
                   "Kohn-Sham fragments and n-mers, and model.method '"// &
                   trim(method_type_to_string(method_type))//"' has no "// &
                   "functional for it to take damping parameters from. Drop "// &
                   "keywords.dft.dispersion, or run model.method 'dft' with a "// &
                   "model.functional."
      end if
   end function dispersion_refusal

   function scheme_name(scheme) result(name)
      !! `scheme` in a sentence: "FMO", "EE-MBE" or "EFMO"
      integer, intent(in) :: scheme
      character(len=:), allocatable :: name

      select case (scheme)
      case (FRAGMENT_SCHEME_FMO)
         name = "FMO"
      case (FRAGMENT_SCHEME_EE_MBE)
         name = "EE-MBE"
      case default
         name = "EFMO"
      end select
   end function scheme_name

   function cut_refusal(scheme, config) result(message)
      !! Why a detached bond is refused for `config%method_type` under `scheme`
      integer, intent(in) :: scheme
      type(method_config_t), intent(in) :: config
      character(len=:), allocatable :: message

      if (scheme == FRAGMENT_SCHEME_EFMO) then
         message = "EFMO: a partition that detaches covalent bonds runs "// &
                   "Hartree-Fock fragments only. The frozen orbitals at a cut are "// &
                   "not yet excluded from the correlation, so an MP2 energy there "// &
                   "would be a different and wrong method; set model.method to 'hf'."
      else
         message = "FMO and EE-MBE: a partition that detaches covalent bonds runs "// &
                   "Hartree-Fock and Kohn-Sham fragments only, and model.method is '"// &
                   trim(method_type_to_string(config%method_type))//"'. The frozen "// &
                   "orbitals at a cut are not yet excluded from the correlation, so a "// &
                   "correlated energy there would be a different and wrong method. Set "// &
                   "keywords.fragmentation.bond_breaking to 'none', or model.method to "// &
                   "'hf' or a functional."
      end if
   end function cut_refusal

end module mqc_fragment_capabilities
