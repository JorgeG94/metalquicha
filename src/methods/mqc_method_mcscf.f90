!! Multi-Configurational Self-Consistent Field (MCSCF) method implementation
module mqc_method_mcscf
   !! Implements CASSCF/CASCI quantum chemistry methods
   !! Provides energy and gradient calculations using complete active space
   !! with optional perturbative corrections (CASPT2/NEVPT2).
   use pic_types, only: dp
   use mqc_method_base, only: qc_method_t
   use mqc_result_types, only: calculation_result_t
   use mqc_physical_fragment, only: physical_fragment_t
   implicit none
   private

   public :: mcscf_method_t, mcscf_options_t

   type :: mcscf_options_t
      !! MCSCF/CASSCF calculation options
      character(len=32) :: basis_set = 'sto-3g'
         !! Basis set name
      logical :: spherical = .true.
         !! Use spherical (true) or Cartesian (false) basis
      logical :: verbose = .false.
         !! Print iteration details

      ! Active space definition
      integer :: n_active_electrons = 0
         !! Number of active electrons (CAS)
      integer :: n_active_orbitals = 0
         !! Number of active orbitals (CAS)
      integer :: n_inactive_orbitals = -1
         !! Number of inactive (doubly occupied) orbitals
         !! -1 means auto-determine from nelec and active electrons

      ! State-averaging
      integer :: n_states = 1
         !! Number of states for state-averaged CASSCF
      real(dp), allocatable :: state_weights(:)
         !! Weights for state averaging (must sum to 1)

      ! Convergence settings
      integer :: max_macro_iter = 100
         !! Maximum macro (orbital optimization) iterations
      integer :: max_micro_iter = 50
         !! Maximum CI iterations per macro step
      real(dp) :: orbital_tol = 1.0e-6_dp
         !! Orbital gradient convergence threshold
      real(dp) :: energy_tol = 1.0e-8_dp
         !! Energy convergence threshold
      real(dp) :: ci_tol = 1.0e-8_dp
         !! CI energy convergence threshold

      ! Orbital optimization algorithm
      character(len=16) :: orbital_optimizer = 'super-ci'
         !! Orbital optimizer: "super-ci", "newton-raphson", "ah" (augmented Hessian)

      ! Perturbative corrections
      logical :: use_pt2 = .false.
         !! Apply perturbative correction after CASSCF
      character(len=16) :: pt2_type = 'nevpt2'
         !! PT2 type: "caspt2", "nevpt2"
      real(dp) :: ipea_shift = 0.25_dp
         !! IPEA shift for CASPT2 (Hartree)
      real(dp) :: imaginary_shift = 0.0_dp
         !! Imaginary shift for intruder states
   end type mcscf_options_t

   type, extends(qc_method_t) :: mcscf_method_t
      !! MCSCF/CASSCF method implementation
      !!
      !! Complete Active Space SCF with optional state-averaging
      !! and perturbative corrections. Suitable for:
      !! - Near-degenerate electronic states
      !! - Bond breaking/formation
      !! - Transition metal complexes
      !! - Excited states
      type(mcscf_options_t) :: options
   contains
      procedure :: calc_energy => mcscf_calc_energy
      procedure :: calc_gradient => mcscf_calc_gradient
      procedure :: calc_hessian => mcscf_calc_hessian
   end type mcscf_method_t

contains

   subroutine mcscf_calc_energy(this, fragment, result)
      !! Calculate electronic energy using CASSCF
      !!
      !! TODO: Implementation requires:
      !! 1. Build basis set and compute integrals
      !! 2. Initial orbital guess (HF or read from file)
      !! 3. Partition orbitals: inactive, active, virtual
      !! 4. Macro iterations:
      !!    a. Transform integrals to MO basis
      !!    b. Solve CI in active space (Davidson or direct)
      !!    c. Build 1- and 2-RDMs from CI vector
      !!    d. Compute orbital gradient
      !!    e. Update orbitals (super-CI, Newton-Raphson, etc.)
      !!    f. Check convergence
      !! 5. Optional: CASPT2/NEVPT2 correction
      class(mcscf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      integer :: n_inactive

      if (this%options%verbose) then
         print *, "MCSCF: Calculating CASSCF energy"
         print *, "MCSCF: Basis set: ", trim(this%options%basis_set)
         print *, "MCSCF: Fragment has", fragment%n_atoms, "atoms"
         print *, "MCSCF: nelec =", fragment%nelec
         print *, "MCSCF: charge =", fragment%charge
         print *, "MCSCF: Active space: (", this%options%n_active_electrons, ",", &
            this%options%n_active_orbitals, ")"

         ! Calculate inactive orbitals
         if (this%options%n_inactive_orbitals < 0) then
            n_inactive = (fragment%nelec - this%options%n_active_electrons)/2
         else
            n_inactive = this%options%n_inactive_orbitals
         end if
         print *, "MCSCF: Inactive orbitals:", n_inactive

         if (this%options%n_states > 1) then
            print *, "MCSCF: State-averaged over", this%options%n_states, "states"
         end if
         if (this%options%use_pt2) then
            print *, "MCSCF: Will apply ", trim(this%options%pt2_type), " correction"
         end if
      end if

      ! Validate active space
      if (this%options%n_active_electrons <= 0 .or. this%options%n_active_orbitals <= 0) then
         print *, "MCSCF: ERROR - Active space not defined!"
         print *, "MCSCF: Set n_active_electrons and n_active_orbitals in config"
         result%has_error = .true.
         return
      end if

      ! Placeholder: Return dummy energy
      ! TODO: Implement actual CASSCF calculation
      result%energy%scf = -1.0_dp*fragment%n_atoms  ! Placeholder
      result%has_energy = .true.

      if (this%options%verbose) then
         print *, "MCSCF: [STUB] CASSCF Energy =", result%energy%total()
      end if

   end subroutine mcscf_calc_energy

   subroutine mcscf_calc_gradient(this, fragment, result)
      !! Calculate energy gradient using CASSCF
      !!
      !! TODO: Implementation requires:
      !! 1. Converged CASSCF (orbitals and CI)
      !! 2. Solve Z-vector (CPHF-like) equations for response
      !! 3. Compute gradient contributions:
      !!    a. One-electron derivative terms
      !!    b. Two-electron derivative terms (with 2-RDM)
      !!    c. Orbital response contribution
      !!    d. CI response contribution (for state-specific)
      !! For state-averaged: gradient of weighted energy
      class(mcscf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      if (this%options%verbose) then
         print *, "MCSCF: Calculating CASSCF gradient"
      end if

      ! First get energy (and converged orbitals/CI)
      call this%calc_energy(fragment, result)
      if (result%has_error) return

      ! Allocate and fill dummy gradient
      allocate (result%gradient(3, fragment%n_atoms))
      result%gradient = 0.0_dp  ! Placeholder
      result%has_gradient = .true.

      if (this%options%verbose) then
         print *, "MCSCF: [STUB] Gradient computed"
      end if

   end subroutine mcscf_calc_gradient

   subroutine mcscf_calc_hessian(this, fragment, result)
      !! Calculate energy Hessian using CASSCF
      !!
      !! TODO: Analytical CASSCF Hessian is very complex:
      !! - Requires second derivatives of integrals
      !! - Coupled-perturbed MCSCF equations
      !! - CI second derivatives
      !! Typically done via finite difference of gradients
      class(mcscf_method_t), intent(in) :: this
      type(physical_fragment_t), intent(in) :: fragment
      type(calculation_result_t), intent(out) :: result

      if (this%options%verbose) then
         print *, "MCSCF: Analytical Hessian not implemented"
         print *, "MCSCF: Use finite difference of gradients instead"
      end if

      ! For now, just compute energy
      call this%calc_energy(fragment, result)
      result%has_hessian = .false.

   end subroutine mcscf_calc_hessian

end module mqc_method_mcscf
