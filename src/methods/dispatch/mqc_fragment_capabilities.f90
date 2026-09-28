!! What a fragmentation scheme may hand a method, and the one refusal site
module mqc_fragment_capabilities
   !! FMO, EE-MBE and EFMO solve their fragments inside the backend rather
   !! than through `qc_method_t`, so what a method may do there is a property
   !! of the pairing of scheme and method. `fragment_capabilities` states it
   !! and `fragment_refusal` is the one place a deck asking for more is turned
   !! away. See `mqc_docs/source/developer_fragment_solver.rst`.
   use mqc_calc_types, only: CALC_TYPE_ENERGY, calc_type_to_string
   use mqc_method_config, only: method_config_t
   use mqc_method_types, only: METHOD_TYPE_HF, METHOD_TYPE_MP2, method_type_to_string
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
   end type fragment_capabilities_t

contains

   pure function fragment_capabilities(scheme, config) result(cap)
      !! What `config%method_type` may do under `scheme`
      !!
      !! FMO and EE-MBE: Hartree-Fock, with or without a detached bond.
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
         end if
      case (FRAGMENT_SCHEME_EFMO)
         select case (config%method_type)
         case (METHOD_TYPE_HF)
            cap%runs = .true.
            cap%cut = .true.
            cap%pieda = .true.
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
      !! all, a non-Energy driver, `unrestricted`, a detached bond, PIEDA.
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
      else if (needs%pieda .and. .not. cap%pieda) then
         message = "keywords.fragmentation.pieda has no decomposition for "// &
                   "model.method '"//trim(method_type_to_string(config%method_type))// &
                   "' yet. Drop pieda, or set model.method to 'hf'."
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
      else
         message = "The fragment calculations of FMO and EE-MBE currently run "// &
                   "Hartree-Fock only, and model.method is '"// &
                   trim(method_type_to_string(config%method_type))//"', which is "// &
                   "not yet wired into them. Set model.method to 'hf'."
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
         message = "FMO and EE-MBE cannot yet detach covalent bonds for "// &
                   "model.method '"//trim(method_type_to_string(config%method_type))// &
                   "'. Set keywords.fragmentation.bond_breaking to 'none', or "// &
                   "model.method to 'hf'."
      end if
   end function cut_refusal

end module mqc_fragment_capabilities
