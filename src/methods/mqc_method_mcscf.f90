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
   !! Everything is computed by `mqc_libcint_mcscf` and `mqc_libcint_casci`
   !! behind `mqc_libcint_bridge`, for the reason `mqc_method_hf` reaches the
   !! CPU SCF the same way: those modules use libcint, fpm globs `src/` and has
   !! no source for it, so `src/methods/` cannot name them. Without that backend
   !! the stub bridge reports the missing build rather than returning a
   !! placeholder -- a plausible wrong number inside a many-body expansion is
   !! far worse than a failed fragment.
   !!
   !! There is no cuEST path and no `#ifdef` here, unlike Hartree-Fock: cuEST
   !! has no CI at all, so there is nothing to choose between.
   use pic_types, only: dp
   use mqc_program_limits, only: MAX_ORBITAL_LABEL_LEN
   use mqc_method_base, only: qc_method_t
   use mqc_method_config, only: pcm_config_t, properties_config_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   use mqc_error, only: ERROR_VALIDATION
   use mqc_cuest_iface, only: apply_properties_settings, cuest_scf_settings_t
   use mqc_libcint_bridge, only: run_libcint_mcscf
   implicit none
   private

   public :: mcscf_method_t, mcscf_options_t

   type :: mcscf_options_t
      !! MCSCF/CASSCF calculation options
      type(properties_config_t) :: properties
      type(pcm_config_t) :: pcm
         !! Continuum solvation. Carried so the refusal can be made where the
         !! request arrives; the CPU path has no cavity.
      character(len=32) :: basis_set = "sto-3g"
         !! Basis set name
      logical :: spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis
      logical :: verbose = .false.
         !! Print iteration details

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
      real(dp) :: energy_tol = 1.0e-8_dp
         !! Energy convergence threshold
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

      ! Reference SCF settings, shared with every other reference-based method
      integer :: max_iter = 100
         !! Maximum reference SCF iterations
      real(dp) :: conv_tol = 1.0e-8_dp
         !! Reference SCF energy convergence threshold
      real(dp) :: density_tol = 1.0e-6_dp
      real(dp) :: linear_dependence = 0.0_dp
         !! Zero means the orthogonaliser's own cutoff. See `scf_config_t`.
         !! Reference SCF density matrix convergence threshold
      real(dp) :: level_shift = 0.0_dp
         !! Hartree added to the virtual block before each diagonalisation.
         !! Zero is off. See `scf_config_t`.
      logical :: use_diis = .true.
      integer :: diis_size = 8
   end type mcscf_options_t

   type, extends(qc_method_t) :: mcscf_method_t
      !! MCSCF/CASSCF method implementation
      !!
      !! Complete Active Space SCF. Suitable for:
      !! - Near-degenerate electronic states
      !! - Bond breaking/formation
      !! - Transition metal complexes
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

      ! The reference SCF's settings, which is most of what this is: a CASSCF
      ! begins with a closed-shell SCF and the backend takes the same settings
      ! object for it that Hartree-Fock does.
      settings%pcm = this%options%pcm

      ! Refused rather than dropped. `run_libcint_mcscf` never reads
      ! `settings%fukui_population` -- unlike `run_libcint_hf`, which does --
      ! so a deck asking for Fukui indices on a CASSCF would otherwise be
      ! obeyed everywhere except where it matters and report a result that
      ! ignored it. That is the failure this whole shared-settings work exists
      ! to end, and handing the field over unread would recreate it one layer
      ! further down.
      if (allocated(this%options%properties%fukui_population)) then
         call result%error%set(ERROR_VALIDATION, "Fukui indices are not available for "// &
                               "an MCSCF reference: they need a response the CASSCF path "// &
                               "does not build. Request them on a Hartree-Fock or "// &
                               "Kohn-Sham calculation instead.")
         result%has_error = .true.
         return
      end if
      call apply_properties_settings(settings, this%options%properties)
      settings%basis_set = this%options%basis_set
      settings%spherical = this%options%spherical
      settings%verbose = this%options%verbose
      settings%functional = ""        ! empty selects pure Hartree-Fock
      settings%max_iter = this%options%max_iter
      settings%energy_tol = this%options%energy_tol
      settings%density_tol = this%options%density_tol
      settings%level_shift = this%options%level_shift
      settings%linear_dependence = this%options%linear_dependence
      settings%use_diis = this%options%use_diis
      settings%diis_size = this%options%diis_size

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

      call run_libcint_mcscf(settings, fragment, result, want_gradient=want_gradient)
   end subroutine mcscf_run

   subroutine mcscf_calc_gradient(this, fragment, result)
      !! The analytic gradient of an optimised CASSCF
      !!
      !! Analytic rather than by differences, and that remains a deliberate
      !! choice: a numerical MCSCF gradient is 6N converged orbital
      !! optimisations, each of which can land on a *different* active space --
      !! the orbitals that make up a CAS are identified by their character, and
      !! a displaced geometry can reorder them. The result has discontinuities
      !! that look like noise and nothing in the output says which displacement
      !! wandered.
      !!
      !! A CASCI is still refused, and the backend says so: its orbitals came
      !! from the SCF and were never optimised for the active space, so the
      !! energy is not stationary with respect to them and the response terms
      !! this omits are not zero.
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
      !! Say that a derivative is not available, once, for both of them
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
