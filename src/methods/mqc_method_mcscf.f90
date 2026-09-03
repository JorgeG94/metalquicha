!! Multi-Configurational Self-Consistent Field (MCSCF) method implementation
module mqc_method_mcscf
   !! CASSCF, and CASCI on the reference orbitals.
   !!
   !! A complete active space wavefunction: the orbitals are partitioned into
   !! doubly occupied *inactive*, an *active* set the CI distributes electrons
   !! over in every possible way, and empty *virtual* ones. CASSCF optimises the
   !! orbitals alongside the CI coefficients; CASCI leaves them at whatever the
   !! reference SCF produced. Both spellings, and "mcscf", parse to
   !! `METHOD_TYPE_MCSCF` -- they are one method with a boolean, not two.
   !!
   !! Everything is computed by `mqc_czt_mcscf` and `mqc_czt_casci`
   !! behind `mqc_czt_bridge`. There is no cuEST path and no `#ifdef` here,
   !! unlike Hartree-Fock: cuEST has no CI at all, so there is nothing to
   !! choose between.
   use pic_types, only: dp
   use mqc_program_limits, only: MAX_ORBITAL_LABEL_LEN
   use mqc_method_base, only: qc_method_t
   use mqc_method_config, only: pcm_config_t, properties_config_t, scf_options_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: ERROR_VALIDATION
   use mqc_cuest_iface, only: apply_properties_settings, apply_scf_settings, &
                              cuest_scf_settings_t
   use mqc_czt_bridge, only: run_czt_mcscf
   implicit none
   private

   public :: mcscf_method_t, mcscf_options_t

   type, extends(scf_options_t) :: mcscf_options_t
      !! What a CASSCF adds to the reference SCF every method shares
      !!
      !! The reference SCF's settings are inherited from `scf_options_t`, not
      !! re-declared. The inherited `energy_tol` is what reaches that SCF.

      ! Active space definition
      character(len=MAX_ORBITAL_LABEL_LEN), allocatable :: avas_orbitals(:)
      real(dp) :: avas_threshold = 0.2_dp
      logical :: full_valence = .false.
      integer :: n_active_electrons = 0
         !! Number of active electrons (CAS)
      integer :: n_active_orbitals = 0
         !! Number of active orbitals (CAS)
      integer :: n_inactive_orbitals = -1
         !! Number of inactive (doubly occupied) orbitals
         !! -1 means auto-determine from nelec and active electrons
      logical :: optimize_orbitals = .true.
         !! CASSCF when true, CASCI when false. See `mcscf_config_t`.
      integer, allocatable :: ormas_subspaces(:)
         !! Active orbital each subspace starts at. Unallocated is a complete
         !! active space.
      integer, allocatable :: ormas_min_electrons(:)
      integer, allocatable :: ormas_max_electrons(:)

      ! State-averaging
      integer :: n_states = 1
         !! Number of states for state-averaged CASSCF
      real(dp), allocatable :: state_weights(:)
         !! Weights for state averaging (must sum to 1)
         !!
         !! Not reachable from a deck: the optimiser underneath solves for one
         !! state, so `mqc_json_schema` allows no key that would set these. The
         !! fields stay because the config type they are copied from has them.

      ! Convergence settings
      integer :: max_macro_iter = 100
         !! Maximum macro (orbital optimization) iterations
      integer :: max_micro_iter = 50
         !! Maximum CI iterations per macro step. Not reachable from a deck:
         !! the Davidson takes a residual threshold, not an iteration cap.
      real(dp) :: orbital_tol = 1.0e-6_dp
         !! Orbital gradient convergence threshold
      real(dp) :: ci_tol = 1.0e-8_dp
         !! CI energy convergence threshold. Not reachable from a deck; the
         !! macro loop pins its inner CASCI tight on purpose.

      ! Orbital optimization algorithm
      character(len=16) :: orbital_optimizer = "super-ci"
         !! Nominal; the implementation is a trust-region Newton step on the
         !! exact orbital Hessian, and there is nothing to select between.

      ! Perturbative corrections
      logical :: use_pt2 = .false.
         !! Apply perturbative correction after CASSCF. No implementation
         !! exists, so no keyword reaches it.
      character(len=16) :: pt2_type = "nevpt2"
         !! PT2 type: "caspt2", "nevpt2"
      real(dp) :: ipea_shift = 0.25_dp
         !! IPEA shift for CASPT2 (Hartree)
      real(dp) :: imaginary_shift = 0.0_dp
         !! Imaginary shift for intruder states

   end type mcscf_options_t

   type, extends(qc_method_t) :: mcscf_method_t
      !! Complete active space SCF: near-degeneracy, bond breaking, open-shell
      !! transition metals
      type(mcscf_options_t) :: options
   contains
      procedure :: calc_energy => mcscf_calc_energy
      procedure :: calc_gradient => mcscf_calc_gradient
      procedure :: calc_hessian => mcscf_calc_hessian
   end type mcscf_method_t

contains

   subroutine mcscf_calc_energy(this, fragment, result)
      !! CASSCF, or CASCI, on one fragment
      class(mcscf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call mcscf_run(this, fragment, result, want_gradient=.false.)
   end subroutine mcscf_calc_energy

   subroutine mcscf_run(this, fragment, result, want_gradient)
      !! The settings the backend needs, built once for both drivers
      class(mcscf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result
      logical, intent(in) :: want_gradient

      type(cuest_scf_settings_t) :: settings

      ! Refused rather than dropped: `run_czt_mcscf` never reads
      ! `settings%fukui_population`, so passing it on would obey the deck
      ! everywhere except where it matters.
      if (allocated(this%options%properties%fukui_population)) then
         call result%error%set(ERROR_VALIDATION, "Fukui indices are not available for "// &
                               "an MCSCF reference: they need a response the CASSCF path "// &
                               "does not build. Request them on a Hartree-Fock or "// &
                               "Kohn-Sham calculation instead.")
         result%has_error = .true.
         return
      end if
      call apply_properties_settings(settings, this%options%properties)
      ! The reference SCF's settings, `pcm` among them: a CASSCF begins with a
      ! closed-shell SCF and the backend takes the same settings object for it.
      call apply_scf_settings(settings, this%options)
      settings%functional = ""        ! empty selects pure Hartree-Fock

      if (allocated(this%options%avas_orbitals)) then
         settings%mcscf%avas_orbitals = this%options%avas_orbitals
      end if
      settings%mcscf%avas_threshold = this%options%avas_threshold
      settings%mcscf%full_valence = this%options%full_valence
      settings%mcscf%n_active_electrons = this%options%n_active_electrons
      settings%mcscf%n_active_orbitals = this%options%n_active_orbitals
      settings%mcscf%n_inactive_orbitals = this%options%n_inactive_orbitals
      settings%mcscf%optimize_orbitals = this%options%optimize_orbitals
      if (allocated(this%options%ormas_subspaces)) then
         settings%mcscf%ormas_subspaces = this%options%ormas_subspaces
         settings%mcscf%ormas_min_electrons = this%options%ormas_min_electrons
         settings%mcscf%ormas_max_electrons = this%options%ormas_max_electrons
      end if
      settings%mcscf%max_macro_iter = this%options%max_macro_iter
      settings%mcscf%orbital_convergence = this%options%orbital_tol

      call run_czt_mcscf(settings, fragment, result, want_gradient=want_gradient)
   end subroutine mcscf_run

   subroutine mcscf_calc_gradient(this, fragment, result)
      !! The analytic gradient of an optimised CASSCF
      !!
      !! Analytic rather than by differences: a displaced geometry can reorder
      !! the orbitals a CAS is built from, so 6N re-optimisations can each land
      !! on a different active space and the result has discontinuities that
      !! look like noise.
      !!
      !! A CASCI is refused by the backend. Its orbitals were never optimised
      !! for the active space, so the energy is not stationary with respect to
      !! them and the response terms this omits are not zero.
      class(mcscf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call mcscf_run(this, fragment, result, want_gradient=.true.)
   end subroutine mcscf_calc_gradient

   subroutine mcscf_calc_hessian(this, fragment, result)
      !! Refused, for the reason the gradient is: there is nothing to difference
      class(mcscf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      call refuse_derivative(fragment, result, "Hessian")
      if (this%options%n_active_orbitals < 0) return
   end subroutine mcscf_calc_hessian

   subroutine refuse_derivative(fragment, result, what)
      !! Say that a derivative is not available
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result
      character(len=*), intent(in) :: what

      call result%error%set(ERROR_VALIDATION, "a CASSCF "//what//" is not implemented: "// &
                            "the analytic form needs the orbital and CI Lagrangian, and "// &
                            "differencing the energy numerically is unreliable here "// &
                            "because the active space can reorder between displacements. "// &
                            "Run driver 'Energy'.")
      result%has_error = .true.
      result%has_energy = .false.
      if (fragment%n_atoms < 0) return
   end subroutine refuse_derivative

end module mqc_method_mcscf
