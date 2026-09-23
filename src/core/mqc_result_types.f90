!! Quantum chemistry calculation result containers
module mqc_result_types
   !! Energy, gradient and property containers, and the MPI transfer of one.
   use pic_types, only: dp, int32, int64
   use pic_mpi_lib, only: comm_t, isend, irecv, send, recv, wait, request_t, MPI_Status
   use mqc_error, only: error_t
   use mqc_calculation_defaults, only: STATE_SPIN_UNKNOWN, STATE_SPIN_SINGLET, &
                                       STATE_SPIN_TRIPLET, STATE_SPIN_UNRESTRICTED
   implicit none
   private

   public :: mp2_energy_t          !! MP2 energy components type
   public :: cc_energy_t           !! Coupled cluster energy components type
   public :: energy_t              !! Energy components type
   public :: calculation_result_t  !! Main result container type
   public :: SCF_UNKNOWN, SCF_CONVERGED, SCF_NOT_CONVERGED
   public :: STATE_SPIN_UNKNOWN, STATE_SPIN_SINGLET, STATE_SPIN_TRIPLET
   public :: STATE_SPIN_UNRESTRICTED
   public :: frontier_orbitals
   public :: scf_not_converged_message
   public :: scf_status_label

   integer, parameter :: SCF_UNKNOWN = 0
      !! The method did not say. Not a claim that it converged.
   integer, parameter :: SCF_CONVERGED = 1
      !! Reached its energy and density thresholds
   integer, parameter :: SCF_NOT_CONVERGED = 2
      !! Ran out of iterations. The energy is whatever the last cycle held --
      !! a number of the right magnitude, so an expansion that sums it is
      !! silently wrong.
   public :: mbe_result_t          !! MBE aggregated result container type
   public :: result_send, result_isend  !! Send result over MPI
   public :: result_recv, result_irecv  !! Receive result over MPI

   ! `STATE_SPIN_*` are defined in `mqc_calculation_defaults` and re-exported
   ! above, `STATE_SPIN_UNRESTRICTED` among them. They are declared beside the
   ! excitation solver's other constants so that a backend can write the label
   ! without importing this module, which carries `pic_mpi_lib` for the result
   ! transfer below.

   ! SCS-MP2 scaling parameters
   real(dp), parameter :: SCS_SS_SCALE = 1.0_dp/3.0_dp  !! SCS same-spin scaling factor
   real(dp), parameter :: SCS_OS_SCALE = 1.2_dp         !! SCS opposite-spin scaling factor

   type :: mp2_energy_t
      !! Container for MP2 energy components (SS/OS)
      real(dp) :: ss = 0.0_dp  !! Same-spin correlation energy (Hartree)
      real(dp) :: os = 0.0_dp  !! Opposite-spin correlation energy (Hartree)
      ! Spin-component scaling, applied by `total` rather than folded into ss
      ! and os, so those stay the numbers the theory produced. Both default to
      ! one, which is plain MP2.
      real(dp) :: ss_scale = 1.0_dp
      real(dp) :: os_scale = 1.0_dp
   contains
      procedure :: total => mp2_total           !! Compute total MP2 correlation
      procedure :: scs => mp2_scs               !! Compute SCS-MP2 correlation
      procedure :: reset => mp2_reset           !! Reset both components to zero
      procedure :: check_stability => mp2_check_stability  !! Check for positive energies (instability)
   end type mp2_energy_t

   type :: cc_energy_t
      !! Container for coupled cluster energy components
      real(dp) :: singles = 0.0_dp   !! Singles contribution (Hartree)
      real(dp) :: doubles = 0.0_dp   !! Doubles contribution (Hartree)
      real(dp) :: triples = 0.0_dp   !! Triples contribution (Hartree)
   contains
      procedure :: total => cc_total            !! Compute total CC correlation
      procedure :: reset => cc_reset            !! Reset all components to zero
      procedure :: check_stability => cc_check_stability  !! Check for positive energies (instability)
   end type cc_energy_t

   type :: energy_t
      !! Container for quantum chemistry energy components
      !!
      !! The total is `scf + dispersion + dh_pt2 + mp2%total() + cc%total()`.
      real(dp) :: scf = 0.0_dp           !! SCF/HF/KS reference energy (Hartree)
      real(dp) :: dispersion = 0.0_dp
         !! Empirical dispersion, from `keywords.dft.dispersion`. Kept out of
         !! `scf` rather than folded into it: it is not a functional of the
         !! density, it converges nothing, and a total that hides it cannot be
         !! compared with a published DFT-D number or with the same geometry run
         !! without it.
      real(dp) :: dh_pt2 = 0.0_dp
         !! The perturbative part of a double hybrid's *functional*.
      type(mp2_energy_t) :: mp2          !! MP2 correlation components
      type(cc_energy_t) :: cc            !! Coupled cluster correlation components
      ! add more as needed, also need to modify the total energy function
   contains
      procedure :: total => energy_total     !! Compute total energy from components
      procedure :: reset => energy_reset     !! Reset all components to zero
   end type energy_t

   type :: calculation_result_t
      !! Container for quantum chemistry calculation results
      !!
      !! Each quantity carries a `has_*` flag saying whether it was filled in.
      type(energy_t) :: energy                  !! Energy components (Hartree)
      real(dp), allocatable :: gradient(:, :)   !! Energy gradient (3, natoms) (Hartree/Bohr)
      real(dp), allocatable :: sigma(:, :)     !! Stress tensor (3,3) (Hartree/Bohr^3)
      real(dp), allocatable :: hessian(:, :)    !! Energy hessian (future implementation)
      real(dp), allocatable :: dipole(:)        !! Dipole moment vector (3) (Debye)
      real(dp), allocatable :: dipole_derivatives(:, :)  !! Dipole derivatives (3, 3N) in a.u. for IR intensities
      real(dp), allocatable :: bond_orders(:, :)
         !! Bond orders, (natoms, natoms), symmetric with a zero diagonal.
         !!
         !! Filled either by an xTB single point (Wiberg-Mayer, in that
         !! Hamiltonian's own minimal basis) or from a converged ab initio
         !! density (Mayer, in the AO basis). Those are different quantities
         !! and only agree to a trend, so `bond_order_scheme` says which one is
         !! here and a consumer that ignores it is comparing two numbers that
         !! were never the same number.
      real(dp), allocatable :: bond_order_valences(:)
         !! (natoms) `sum_B B_AB`, the valence each atom's orders add up to.
         !! Free once the matrix exists and the column a reader looks at first.
      character(len=16) :: bond_order_scheme = ""
         !! Which definition produced them: "mayer" or "xtb". Empty means
         !! nothing has.

      ! Atomic partial charges, when `properties.charges` asked for them.
      real(dp), allocatable :: atomic_charges(:)
         !! (natoms) nuclear charge minus the electrons assigned to the atom, so
         !! positive is electron-poor.
      real(dp), allocatable :: spin_populations(:)
         !! (natoms) the same partition applied to `P_alpha - P_beta`, so it
         !! sums to `n_alpha - n_beta`. Allocated only for an unrestricted
         !! reference under Mulliken; CHELPG fits the electrostatic potential,
         !! which has no spin analogue.
      character(len=16) :: charge_scheme = ""
         !! Which partition produced them. Two schemes disagree by design, so a
         !! charge without its scheme is not interpretable.
      logical :: has_charges = .false.

      ! Intrinsic energy decomposition, when `properties.bonding_analysis`
      ! asked for one. Hartree throughout, as everything here is; the printed
      ! tables convert.

      real(dp), allocatable :: ieda_atom(:)
         !! (natoms) everything that happens inside one atom -- its own kinetic
         !! energy, its electrons in its own nuclear field, its own repulsion.
         !! An atom deformed by the molecule, not a free one.
      real(dp), allocatable :: ieda_free_atom(:)
         !! (natoms) the same atom on its own. `ieda_atom - ieda_free_atom` is
         !! what it cost to prepare it in the shape the molecule needs.
      real(dp), allocatable :: ieda_pair(:, :)
         !! (natoms, natoms) everything that exists only because two atoms are
         !! both present. **The full pair energy sits in both (A,B) and (B,A)**,
         !! which is the convention the decomposition itself uses, so the total
         !! is `sum(ieda_atom) + 0.5*sum(ieda_pair)` and not a triangle sum.
      real(dp), allocatable :: ieda_classical(:, :)
         !! (natoms, natoms) the part of the pair energy an electrostatic model
         !! could produce; the remainder is interference.
      real(dp) :: ieda_formation = 0.0_dp
         !! The energy of formation: the molecule against its free atoms.
      logical :: has_ieda = .false.
      real(dp), allocatable :: fukui_plus(:)
         !! Condensed Fukui index for nucleophilic attack, per atom
      real(dp), allocatable :: fukui_minus(:)   !! ... for electrophilic attack
      real(dp), allocatable :: fukui_dual(:)    !! `f+ - f-`
      real(dp) :: fukui_ip = 0.0_dp             !! E(N-1) - E(N)
      real(dp) :: fukui_ea = 0.0_dp             !! E(N) - E(N+1)
      real(dp) :: fukui_hardness = 0.0_dp
      real(dp) :: fukui_electrophilicity = 0.0_dp
      logical :: fukui_anion_bound = .true.
         !! False when the anion is unbound, in which case `f+` describes an
         !! orbital the basis invented and nothing about the numbers says so.
      character(len=16) :: fukui_scheme = ""
      logical :: has_fukui = .false.

      ! Linear-response excited states, when `keywords.excited_states` asked
      ! for any. Every per-state array below runs over the same states in the
      ! same order, lowest excitation first. The energies and spins arrive
      ! with the spectrum; the moments arrive with the properties, which are
      ! computed after it and can fail on their own, so a spectrum without
      ! them is a state this type has to be able to hold.
      real(dp), allocatable :: excitation_energies(:)
         !! (n_states) vertical excitation energies in Hartree, above the
         !! reference total energy rather than absolute.
      real(dp), allocatable :: excited_total_energies(:)
         !! (n_states) each excited state's own total energy, `E_SCF + w`, in
         !! Hartree. Carried rather than left to the consumer because the
         !! reference energy a fragment's spectrum sits on is that fragment's,
         !! not the one printed beside it.
      real(dp), allocatable :: oscillator_strengths(:)
         !! (n_states) dimensionless length-gauge oscillator strengths. Exactly
         !! zero for a triplet, which is a real value and not a missing one.
      real(dp), allocatable :: oscillator_strengths_velocity(:)
         !! (n_states) the same strengths in the velocity gauge. The two
         !! agree only in a complete basis, so both are reported: the gap
         !! between them measures the basis rather than the solve.
      real(dp), allocatable :: transition_dipoles(:, :)
         !! (3, n_states) transition dipole moments in atomic units, with the
         !! origin at the nuclear charge centroid.
      real(dp), allocatable :: transition_velocities(:, :)
         !! (3, n_states) the same transitions in the velocity gauge, atomic
         !! units: the imaginary part of `<0|p|n>`. The length-gauge partner
         !! of a row is the one directly above it, and comparing the two is
         !! what the second gauge is reported for.
      real(dp), allocatable :: transition_dipole_origin(:)
         !! (3) where the length-gauge dipoles were measured from, in Bohr.
         !! One origin for the whole spectrum. A transition dipole is origin
         !! independent for a neutral transition density, so this says which
         !! convention was used rather than changing the numbers.
      real(dp), allocatable :: nto_leading_weight(:)
         !! (n_states) the largest natural transition orbital weight, between
         !! zero and one. One says the root is exactly one orbital pair.
      integer, allocatable :: state_spin(:)
         !! (n_states) which spin each root carries, as `STATE_SPIN_*`. Carried
         !! per state rather than once for the run because a `spin: "both"`
         !! deck interleaves singlets and triplets in one ordered list, and
         !! nothing in an excitation energy says which kind it is.
      logical :: has_excited_states = .false.

      logical :: stability_stable = .true.
         !! Whether the converged SCF is a minimum with respect to real
         !! closed-shell orbital rotations. Meaningful only with
         !! `has_stability`; see `mqc_czt_ov_hessian` for what it does not
         !! cover -- a triplet or complex instability is a different matrix.
      real(dp) :: stability_curvature = 0.0_dp
         !! The lowest eigenvalue of the electronic Hessian, in hartree.
      logical :: stability_has_curvature = .false.
         !! Whether that eigenvalue was recovered at all. It always is on
         !! this path, and the flag is in the output contract anyway so that a
         !! consumer can tell an absent number from an absent feature.
      integer :: stability_rotations = 0
         !! How many non-redundant orbital rotations were searched.
      logical :: has_stability = .false.

      ! Fragment metadata
      real(dp) :: distance = 0.0_dp      !! Minimal atomic distance between monomers (Angstrom, 0 for monomers)

      ! Computation status flags
      logical :: has_energy = .false.    !! Energy has been computed
      logical :: has_gradient = .false.  !! Gradient has been computed
      logical :: has_sigma = .false.     !! Stress tensor has been computed
      logical :: has_hessian = .false.   !! Hessian has been computed
      logical :: has_dipole = .false.    !! Dipole moment has been computed
      logical :: has_dipole_derivatives = .false.  !! Dipole derivatives have been computed
      logical :: has_bond_orders = .false.  !! Bond orders have been computed

      ! Frontier orbitals only, not the whole spectrum
      real(dp) :: homo = 0.0_dp          !! Highest occupied orbital energy (Hartree)
      real(dp) :: lumo = 0.0_dp          !! Lowest unoccupied orbital energy (Hartree)
      logical :: has_orbitals = .false.  !! Whether homo/lumo were reported

      ! SCF convergence. Tri-state rather than a logical: a method that does
      ! not report leaves this at SCF_UNKNOWN and the output prints "?",
      ! neither claiming convergence nor claiming failure.
      integer :: scf_status = SCF_UNKNOWN  !! Whether the SCF reached its thresholds
      integer :: scf_iterations = 0        !! Cycles used, 0 if not reported

      ! Error handling
      type(error_t) :: error             !! Calculation error (if any)
      logical :: has_error = .false.     !! True if calculation failed
   contains
      procedure :: destroy => result_destroy  !! Clean up allocated memory
      procedure :: reset => result_reset      !! Reset all values and flags
   end type calculation_result_t

   type :: mbe_result_t
      !! Container for Many-Body Expansion aggregated results
      !!
      !! The caller allocates the components it wants before calling
      !! `compute_mbe`, which computes whatever is allocated and sets the
      !! matching `has_*` flag on success.

      real(dp) :: total_energy = 0.0_dp              !! Total MBE energy (Hartree)
      real(dp), allocatable :: gradient(:, :)        !! Total gradient (3, total_atoms) (Hartree/Bohr)
      real(dp), allocatable :: hessian(:, :)         !! Total Hessian (3*natoms, 3*natoms)
      real(dp), allocatable :: dipole(:)             !! Total dipole moment (3) (e*Bohr)
      real(dp), allocatable :: dipole_derivatives(:, :)  !! Dipole derivatives (3, 3*natoms) for IR intensities

      ! Computation status flags
      logical :: has_energy = .false.                !! Energy has been computed
      logical :: has_gradient = .false.              !! Gradient has been computed
      logical :: has_hessian = .false.               !! Hessian has been computed
      logical :: has_dipole = .false.                !! Dipole has been computed
      logical :: has_dipole_derivatives = .false.    !! Dipole derivatives have been computed

      ! The interaction energy of one fragment, when `compute_mbe` was given
      ! a reference. `total_energy` and `has_energy` are then left unset: the
      ! expansion was reduced to the terms these need, and what it sums to is
      ! not the system's energy.
      integer :: reference_fragment = 0
         !! The reference as a monomer number, 1-based; 0 when there was none
      real(dp) :: reference_energy = 0.0_dp
         !! The reference fragment's own energy, its one-body term (Hartree)
      real(dp) :: interaction_energy = 0.0_dp
         !! Sum of every many-body correction, level 2 and up, whose term
         !! contains the reference (Hartree)
      real(dp), allocatable :: interaction_by_level(:)
         !! (max_level) the same sum split by term size; element 1 is zero
      integer(int64), allocatable :: interaction_count_by_level(:)
         !! (max_level) how many terms containing the reference each level has
      logical :: has_interaction = .false.           !! The four above are set
   contains
      procedure :: destroy => mbe_result_destroy            !! Clean up allocated memory
      procedure :: reset => mbe_result_reset                !! Reset all values and flags
      procedure :: allocate_gradient => mbe_result_allocate_gradient
      procedure :: allocate_hessian => mbe_result_allocate_hessian
      procedure :: allocate_dipole => mbe_result_allocate_dipole
   end type mbe_result_t

contains

   pure function mp2_total(this) result(total)
      !! The correlation energy as scaled, which for plain MP2 is the sum
      class(mp2_energy_t), intent(in) :: this
      real(dp) :: total

      total = this%ss_scale*this%ss + this%os_scale*this%os
   end function mp2_total

   pure function mp2_scs(this) result(scs_energy)
      !! Grimme's SCS-MP2 at its published factors, whatever this run used
      !!
      !! The factors are fixed here, where `total` applies the run's own.
      class(mp2_energy_t), intent(in) :: this
      real(dp) :: scs_energy

      scs_energy = SCS_SS_SCALE*this%ss + SCS_OS_SCALE*this%os
   end function mp2_scs

   subroutine mp2_reset(this)
      !! Reset both MP2 components to zero
      class(mp2_energy_t), intent(inout) :: this
      this%ss = 0.0_dp
      this%os = 0.0_dp
   end subroutine mp2_reset

   subroutine mp2_check_stability(this)
      !! Check for positive MP2 correlation energies (instability warning)
      !! Correlation energies should be negative; positive values indicate instability
      use pic_logger, only: logger => global_logger
      class(mp2_energy_t), intent(in) :: this

      if (this%ss > 0.0_dp) then
         call logger%warning("MP2 same-spin correlation energy is positive - possible instability!")
      end if

      if (this%os > 0.0_dp) then
         call logger%warning("MP2 opposite-spin correlation energy is positive - possible instability!")
      end if
   end subroutine mp2_check_stability

   pure function cc_total(this) result(total)
      !! Compute total CC correlation energy
      class(cc_energy_t), intent(in) :: this
      real(dp) :: total

      total = this%singles + this%doubles + this%triples
   end function cc_total

   subroutine cc_reset(this)
      !! Reset all CC components to zero
      class(cc_energy_t), intent(inout) :: this
      this%singles = 0.0_dp
      this%doubles = 0.0_dp
      this%triples = 0.0_dp
   end subroutine cc_reset

   subroutine cc_check_stability(this)
      !! Check for positive CC correlation energies (instability warning)
      !! Correlation energies should be negative; positive values indicate instability
      use pic_logger, only: logger => global_logger
      class(cc_energy_t), intent(in) :: this

      if (this%singles > 0.0_dp) then
         call logger%warning("CC singles correlation energy is positive - possible instability!")
      end if

      if (this%doubles > 0.0_dp) then
         call logger%warning("CC doubles correlation energy is positive - possible instability!")
      end if

      if (this%triples > 0.0_dp) then
         call logger%warning("CC triples correlation energy is positive - possible instability!")
      end if
   end subroutine cc_check_stability

   pure function energy_total(this) result(total)
      !! Compute total energy from all components
      class(energy_t), intent(in) :: this
      real(dp) :: total

      ! this line needs to me modified if more components are added
      total = this%scf + this%dispersion + this%dh_pt2 + this%mp2%total() + this%cc%total()
   end function energy_total

   subroutine energy_reset(this)
      !! Reset all energy components to zero
      class(energy_t), intent(inout) :: this
      this%scf = 0.0_dp
      this%dispersion = 0.0_dp
      this%dh_pt2 = 0.0_dp
      call this%mp2%reset()
      call this%cc%reset()
   end subroutine energy_reset

   subroutine result_destroy(this)
      !! Clean up allocated memory in calculation_result_t
      ! TODO(mqc): `bond_orders`, `fukui_plus`, `fukui_minus` and `fukui_dual`
      ! are never deallocated here, so a result reused through `result_recv`
      ! keeps the previous fragment's arrays at the previous fragment's size
      ! while `has_bond_orders` still says they are current.
      class(calculation_result_t), intent(inout) :: this
      if (allocated(this%gradient)) deallocate (this%gradient)
      if (allocated(this%sigma)) deallocate (this%sigma)
      if (allocated(this%hessian)) deallocate (this%hessian)
      if (allocated(this%dipole)) deallocate (this%dipole)
      if (allocated(this%dipole_derivatives)) deallocate (this%dipole_derivatives)
      if (allocated(this%atomic_charges)) deallocate (this%atomic_charges)
      if (allocated(this%spin_populations)) deallocate (this%spin_populations)
      if (allocated(this%ieda_atom)) deallocate (this%ieda_atom)
      if (allocated(this%ieda_free_atom)) deallocate (this%ieda_free_atom)
      if (allocated(this%ieda_pair)) deallocate (this%ieda_pair)
      if (allocated(this%ieda_classical)) deallocate (this%ieda_classical)
      ! `result_recv` destroys precisely so a container can be reused, so an
      ! array left behind here comes back on the next fragment at the previous
      ! fragment's size.
      if (allocated(this%bond_orders)) deallocate (this%bond_orders)
      if (allocated(this%bond_order_valences)) deallocate (this%bond_order_valences)
      if (allocated(this%fukui_plus)) deallocate (this%fukui_plus)
      if (allocated(this%fukui_minus)) deallocate (this%fukui_minus)
      if (allocated(this%fukui_dual)) deallocate (this%fukui_dual)
      if (allocated(this%excitation_energies)) deallocate (this%excitation_energies)
      if (allocated(this%excited_total_energies)) then
         deallocate (this%excited_total_energies)
      end if
      if (allocated(this%oscillator_strengths)) deallocate (this%oscillator_strengths)
      if (allocated(this%oscillator_strengths_velocity)) then
         deallocate (this%oscillator_strengths_velocity)
      end if
      if (allocated(this%transition_dipoles)) deallocate (this%transition_dipoles)
      if (allocated(this%transition_velocities)) then
         deallocate (this%transition_velocities)
      end if
      if (allocated(this%transition_dipole_origin)) then
         deallocate (this%transition_dipole_origin)
      end if
      if (allocated(this%nto_leading_weight)) deallocate (this%nto_leading_weight)
      if (allocated(this%state_spin)) deallocate (this%state_spin)
      call this%reset()
   end subroutine result_destroy

   subroutine result_reset(this)
      !! Reset all values and flags in calculation_result_t
      ! TODO(mqc): clears neither `has_charges` nor `has_bond_orders`, and
      ! `destroy` deallocates `atomic_charges` on the way through, so a reused
      ! result claims charges it no longer holds.
      class(calculation_result_t), intent(inout) :: this
      call this%energy%reset()
      call this%error%clear()
      this%has_energy = .false.
      this%has_gradient = .false.
      this%has_sigma = .false.
      this%has_hessian = .false.
      this%has_dipole = .false.
      this%has_dipole_derivatives = .false.
      this%has_error = .false.
      this%scf_status = SCF_UNKNOWN
      this%scf_iterations = 0
      this%homo = 0.0_dp
      this%lumo = 0.0_dp
      this%has_orbitals = .false.
      this%has_ieda = .false.
      this%has_fukui = .false.
      this%has_excited_states = .false.
      this%has_stability = .false.
      this%stability_stable = .true.
      this%stability_has_curvature = .false.
      this%stability_curvature = 0.0_dp
      this%stability_rotations = 0
      this%ieda_formation = 0.0_dp
   end subroutine result_reset

   !---------------------------------------------------------------------------
   ! mbe_result_t type-bound procedures
   !---------------------------------------------------------------------------

   subroutine mbe_result_destroy(this)
      !! Clean up allocated memory in mbe_result_t
      class(mbe_result_t), intent(inout) :: this
      if (allocated(this%gradient)) deallocate (this%gradient)
      if (allocated(this%hessian)) deallocate (this%hessian)
      if (allocated(this%dipole)) deallocate (this%dipole)
      if (allocated(this%dipole_derivatives)) deallocate (this%dipole_derivatives)
      if (allocated(this%interaction_by_level)) deallocate (this%interaction_by_level)
      if (allocated(this%interaction_count_by_level)) deallocate (this%interaction_count_by_level)
      call this%reset()
   end subroutine mbe_result_destroy

   subroutine mbe_result_reset(this)
      !! Reset all values and flags in mbe_result_t
      class(mbe_result_t), intent(inout) :: this
      this%total_energy = 0.0_dp
      this%has_energy = .false.
      this%has_gradient = .false.
      this%has_hessian = .false.
      this%has_dipole = .false.
      this%has_dipole_derivatives = .false.
      this%reference_fragment = 0
      this%reference_energy = 0.0_dp
      this%interaction_energy = 0.0_dp
      this%has_interaction = .false.
   end subroutine mbe_result_reset

   subroutine mbe_result_allocate_gradient(this, total_atoms)
      !! Allocate gradient array for total_atoms
      class(mbe_result_t), intent(inout) :: this
      integer, intent(in) :: total_atoms
      if (allocated(this%gradient)) deallocate (this%gradient)
      allocate (this%gradient(3, total_atoms))
      this%gradient = 0.0_dp
   end subroutine mbe_result_allocate_gradient

   subroutine mbe_result_allocate_hessian(this, total_atoms)
      !! Allocate hessian array for total_atoms
      class(mbe_result_t), intent(inout) :: this
      integer, intent(in) :: total_atoms
      integer :: hess_dim
      hess_dim = 3*total_atoms
      if (allocated(this%hessian)) deallocate (this%hessian)
      allocate (this%hessian(hess_dim, hess_dim))
      this%hessian = 0.0_dp
   end subroutine mbe_result_allocate_hessian

   subroutine mbe_result_allocate_dipole(this)
      !! Allocate dipole array (always 3 components)
      class(mbe_result_t), intent(inout) :: this
      if (allocated(this%dipole)) deallocate (this%dipole)
      allocate (this%dipole(3))
      this%dipole = 0.0_dp
   end subroutine mbe_result_allocate_dipole

   subroutine send_error_state(result, comm, dest, tag)
      !! Send the failure flag, code and trace that travel with a result
      !!
      !! Every `has_error` check on the receiving side depends on this being
      !! sent. Without it a fragment that failed on a worker arrives looking
      !! healthy and its untouched SCF energy of 0 is summed at face value.
      type(calculation_result_t), intent(in) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: dest, tag

      character(len=:), allocatable :: trace
      integer(int32), allocatable :: encoded(:)
      integer :: i, n

      call send(comm, result%has_error, dest, tag)
      if (.not. result%has_error) return

      call send(comm, result%error%get_code(), dest, tag)

      ! pic-mpi carries no character type, so the trace goes as one character
      ! code per element, and is never empty -- a zero-length message would
      ! leave the coordinator with a failure it cannot name.
      trace = result%error%get_full_trace()
      if (len_trim(trace) == 0) trace = "Calculation failed without a diagnostic"
      n = len_trim(trace)
      allocate (encoded(n))
      do i = 1, n
         encoded(i) = iachar(trace(i:i))
      end do
      call send(comm, encoded, dest, tag)
   end subroutine send_error_state

   subroutine recv_error_state(result, comm, source, tag)
      !! Receive the failure flag, code and trace that travel with a result
      type(calculation_result_t), intent(inout) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: source, tag

      type(MPI_Status) :: status
      integer(int32), allocatable :: encoded(:)
      character(len=:), allocatable :: trace
      integer(int32) :: code
      integer :: i

      call recv(comm, result%has_error, source, tag, status)
      if (.not. result%has_error) then
         call result%error%clear()
         return
      end if

      call recv(comm, code, source, tag, status)
      call recv(comm, encoded, source, tag, status)

      allocate (character(len=size(encoded)) :: trace)
      do i = 1, size(encoded)
         trace(i:i) = achar(encoded(i))
      end do
      call result%error%set(int(code), trace)

      ! A failed fragment has no usable energy, whatever the flags said.
      result%has_energy = .false.
   end subroutine recv_error_state

   subroutine result_send(result, comm, dest, tag)
      !! Send calculation result over MPI (blocking)
      type(calculation_result_t), intent(in) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: dest, tag

      ! Send energy components
      call send(comm, result%energy%scf, dest, tag)
      call send(comm, result%energy%dispersion, dest, tag)
      call send(comm, result%energy%mp2%ss, dest, tag)
      call send(comm, result%energy%mp2%os, dest, tag)
      call send(comm, result%energy%cc%singles, dest, tag)
      call send(comm, result%energy%cc%doubles, dest, tag)
      call send(comm, result%energy%cc%triples, dest, tag)

      ! Send fragment metadata. A field added to the type and not to these
      ! four routines leaves every parallel run reporting the default.
      call send(comm, result%distance, dest, tag)
      call send(comm, result%scf_status, dest, tag)
      call send(comm, result%scf_iterations, dest, tag)
      call send(comm, result%has_orbitals, dest, tag)
      call send(comm, result%homo, dest, tag)
      call send(comm, result%lumo, dest, tag)

      ! Send gradient flag and data if present
      call send(comm, result%has_gradient, dest, tag)
      if (result%has_gradient) then
         call send(comm, result%gradient, dest, tag)
      end if

      ! Send Hessian flag and data if present
      call send(comm, result%has_hessian, dest, tag)
      if (result%has_hessian) then
         call send(comm, result%hessian, dest, tag)
      end if

      ! Send dipole flag and data if present
      call send(comm, result%has_dipole, dest, tag)
      if (result%has_dipole) then
         call send(comm, result%dipole, dest, tag)
      end if

      ! Send dipole derivatives flag and data if present
      call send(comm, result%has_dipole_derivatives, dest, tag)
      if (result%has_dipole_derivatives) then
         call send(comm, result%dipole_derivatives, dest, tag)
      end if

      ! Excited states, every array under one flag: the solver and the
      ! properties fill them together, so a receiver that got the energies
      ! without the spins could not label a single root.
      call send(comm, result%has_excited_states, dest, tag)
      if (result%has_excited_states) then
         call send(comm, allocated(result%excitation_energies), dest, tag)
         if (allocated(result%excitation_energies)) then
            call send(comm, result%excitation_energies, dest, tag)
         end if
         call send(comm, allocated(result%excited_total_energies), dest, tag)
         if (allocated(result%excited_total_energies)) then
            call send(comm, result%excited_total_energies, dest, tag)
         end if
         call send(comm, allocated(result%oscillator_strengths), dest, tag)
         if (allocated(result%oscillator_strengths)) then
            call send(comm, result%oscillator_strengths, dest, tag)
         end if
         call send(comm, allocated(result%oscillator_strengths_velocity), dest, tag)
         if (allocated(result%oscillator_strengths_velocity)) then
            call send(comm, result%oscillator_strengths_velocity, dest, tag)
         end if
         call send(comm, allocated(result%transition_dipoles), dest, tag)
         if (allocated(result%transition_dipoles)) then
            call send(comm, result%transition_dipoles, dest, tag)
         end if
         call send(comm, allocated(result%transition_velocities), dest, tag)
         if (allocated(result%transition_velocities)) then
            call send(comm, result%transition_velocities, dest, tag)
         end if
         call send(comm, allocated(result%transition_dipole_origin), dest, tag)
         if (allocated(result%transition_dipole_origin)) then
            call send(comm, result%transition_dipole_origin, dest, tag)
         end if
         call send(comm, allocated(result%nto_leading_weight), dest, tag)
         if (allocated(result%nto_leading_weight)) then
            call send(comm, result%nto_leading_weight, dest, tag)
         end if
         call send(comm, allocated(result%state_spin), dest, tag)
         if (allocated(result%state_spin)) then
            call send(comm, result%state_spin, dest, tag)
         end if
      end if

      ! Failure state last, so the receiver has the whole payload drained
      ! before it decides whether to trust any of it.
      call send_error_state(result, comm, dest, tag)
   end subroutine result_send

   subroutine result_isend(result, comm, dest, tag, req)
      !! Send calculation result over MPI (non-blocking)
      !! Sends SCF energy (non-blocking) and other components (blocking)
      type(calculation_result_t), intent(in) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: dest, tag
      type(request_t), intent(out) :: req

      ! Send SCF energy (non-blocking)
      call isend(comm, result%energy%scf, dest, tag, req)

      ! Send other energy components (blocking to avoid needing multiple request handles)
      call send(comm, result%energy%dispersion, dest, tag)
      call send(comm, result%energy%mp2%ss, dest, tag)
      call send(comm, result%energy%mp2%os, dest, tag)
      call send(comm, result%energy%cc%singles, dest, tag)
      call send(comm, result%energy%cc%doubles, dest, tag)
      call send(comm, result%energy%cc%triples, dest, tag)

      ! Send fragment metadata. A field added to the type and not to these
      ! four routines leaves every parallel run reporting the default.
      call send(comm, result%distance, dest, tag)
      call send(comm, result%scf_status, dest, tag)
      call send(comm, result%scf_iterations, dest, tag)
      call send(comm, result%has_orbitals, dest, tag)
      call send(comm, result%homo, dest, tag)
      call send(comm, result%lumo, dest, tag)

      ! Send gradient flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_gradient, dest, tag)
      if (result%has_gradient) then
         call send(comm, result%gradient, dest, tag)
      end if

      ! Send Hessian flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_hessian, dest, tag)
      if (result%has_hessian) then
         call send(comm, result%hessian, dest, tag)
      end if

      ! Send dipole flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_dipole, dest, tag)
      if (result%has_dipole) then
         call send(comm, result%dipole, dest, tag)
      end if

      ! Send dipole derivatives flag and data (blocking to avoid needing multiple request handles)
      call send(comm, result%has_dipole_derivatives, dest, tag)
      if (result%has_dipole_derivatives) then
         call send(comm, result%dipole_derivatives, dest, tag)
      end if

      ! Excited states, every array under one flag: the solver and the
      ! properties fill them together, so a receiver that got the energies
      ! without the spins could not label a single root.
      call send(comm, result%has_excited_states, dest, tag)
      if (result%has_excited_states) then
         call send(comm, allocated(result%excitation_energies), dest, tag)
         if (allocated(result%excitation_energies)) then
            call send(comm, result%excitation_energies, dest, tag)
         end if
         call send(comm, allocated(result%excited_total_energies), dest, tag)
         if (allocated(result%excited_total_energies)) then
            call send(comm, result%excited_total_energies, dest, tag)
         end if
         call send(comm, allocated(result%oscillator_strengths), dest, tag)
         if (allocated(result%oscillator_strengths)) then
            call send(comm, result%oscillator_strengths, dest, tag)
         end if
         call send(comm, allocated(result%oscillator_strengths_velocity), dest, tag)
         if (allocated(result%oscillator_strengths_velocity)) then
            call send(comm, result%oscillator_strengths_velocity, dest, tag)
         end if
         call send(comm, allocated(result%transition_dipoles), dest, tag)
         if (allocated(result%transition_dipoles)) then
            call send(comm, result%transition_dipoles, dest, tag)
         end if
         call send(comm, allocated(result%transition_velocities), dest, tag)
         if (allocated(result%transition_velocities)) then
            call send(comm, result%transition_velocities, dest, tag)
         end if
         call send(comm, allocated(result%transition_dipole_origin), dest, tag)
         if (allocated(result%transition_dipole_origin)) then
            call send(comm, result%transition_dipole_origin, dest, tag)
         end if
         call send(comm, allocated(result%nto_leading_weight), dest, tag)
         if (allocated(result%nto_leading_weight)) then
            call send(comm, result%nto_leading_weight, dest, tag)
         end if
         call send(comm, allocated(result%state_spin), dest, tag)
         if (allocated(result%state_spin)) then
            call send(comm, result%state_spin, dest, tag)
         end if
      end if

      ! Failure state last, so the receiver has the whole payload drained
      ! before it decides whether to trust any of it.
      call send_error_state(result, comm, dest, tag)
   end subroutine result_isend

   subroutine result_recv(result, comm, source, tag, status)
      !! Receive calculation result over MPI (blocking)
      type(calculation_result_t), intent(inout) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: source, tag
      type(MPI_Status), intent(out) :: status
      logical :: present_excitation_energies, present_excited_total_energies
      logical :: present_oscillator_strengths, present_oscillator_strengths_velocity
      logical :: present_transition_dipoles, present_transition_velocities
      logical :: present_transition_dipole_origin, present_nto_leading_weight
      logical :: present_state_spin
         !! Presence bytes for the excited-state arrays, read off the
         !! wire ahead of each one. A spectrum can arrive without its
         !! transition properties; see the note at the send.

      ! Same reuse hazard as result_irecv -- see the note there.
      call result%destroy()

      ! Receive energy components
      call recv(comm, result%energy%scf, source, tag, status)
      call recv(comm, result%energy%dispersion, source, tag, status)
      call recv(comm, result%energy%mp2%ss, source, tag, status)
      call recv(comm, result%energy%mp2%os, source, tag, status)
      call recv(comm, result%energy%cc%singles, source, tag, status)
      call recv(comm, result%energy%cc%doubles, source, tag, status)
      call recv(comm, result%energy%cc%triples, source, tag, status)
      result%has_energy = .true.

      ! Receive fragment metadata
      call recv(comm, result%distance, source, tag, status)
      call recv(comm, result%scf_status, source, tag, status)
      call recv(comm, result%scf_iterations, source, tag, status)
      call recv(comm, result%has_orbitals, source, tag, status)
      call recv(comm, result%homo, source, tag, status)
      call recv(comm, result%lumo, source, tag, status)

      ! Receive gradient flag and data if present
      call recv(comm, result%has_gradient, source, tag, status)
      if (result%has_gradient) then
         ! Receive allocatable gradient array (MPI lib handles allocation)
         call recv(comm, result%gradient, source, tag, status)
      end if

      ! Receive Hessian flag and data if present
      call recv(comm, result%has_hessian, source, tag, status)
      if (result%has_hessian) then
         ! Receive allocatable Hessian array (MPI lib handles allocation)
         call recv(comm, result%hessian, source, tag, status)
      end if

      ! Receive dipole flag and data if present
      call recv(comm, result%has_dipole, source, tag, status)
      if (result%has_dipole) then
         ! Receive allocatable dipole array (MPI lib handles allocation)
         call recv(comm, result%dipole, source, tag, status)
      end if

      ! Receive dipole derivatives flag and data if present
      call recv(comm, result%has_dipole_derivatives, source, tag, status)
      if (result%has_dipole_derivatives) then
         ! Receive allocatable dipole derivatives array (MPI lib handles allocation)
         call recv(comm, result%dipole_derivatives, source, tag, status)
      end if

      ! Receive excited states, every array under one flag
      call recv(comm, result%has_excited_states, source, tag, status)
      if (result%has_excited_states) then
         call recv(comm, present_excitation_energies, source, tag, status)
         if (present_excitation_energies) then
            call recv(comm, result%excitation_energies, source, tag, status)
         end if
         call recv(comm, present_excited_total_energies, source, tag, status)
         if (present_excited_total_energies) then
            call recv(comm, result%excited_total_energies, source, tag, status)
         end if
         call recv(comm, present_oscillator_strengths, source, tag, status)
         if (present_oscillator_strengths) then
            call recv(comm, result%oscillator_strengths, source, tag, status)
         end if
         call recv(comm, present_oscillator_strengths_velocity, source, tag, status)
         if (present_oscillator_strengths_velocity) then
            call recv(comm, result%oscillator_strengths_velocity, source, tag, status)
         end if
         call recv(comm, present_transition_dipoles, source, tag, status)
         if (present_transition_dipoles) then
            call recv(comm, result%transition_dipoles, source, tag, status)
         end if
         call recv(comm, present_transition_velocities, source, tag, status)
         if (present_transition_velocities) then
            call recv(comm, result%transition_velocities, source, tag, status)
         end if
         call recv(comm, present_transition_dipole_origin, source, tag, status)
         if (present_transition_dipole_origin) then
            call recv(comm, result%transition_dipole_origin, source, tag, status)
         end if
         call recv(comm, present_nto_leading_weight, source, tag, status)
         if (present_nto_leading_weight) then
            call recv(comm, result%nto_leading_weight, source, tag, status)
         end if
         call recv(comm, present_state_spin, source, tag, status)
         if (present_state_spin) then
            call recv(comm, result%state_spin, source, tag, status)
         end if
      end if

      call recv_error_state(result, comm, source, tag)
   end subroutine result_recv

   subroutine result_irecv(result, comm, source, tag, req)
      !! Receive calculation result over MPI (non-blocking)
      !! Receives SCF energy (non-blocking) and other components (blocking)
      type(calculation_result_t), intent(inout) :: result
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: source, tag
      type(request_t), intent(out) :: req
      type(MPI_Status) :: status
      logical :: present_excitation_energies, present_excited_total_energies
      logical :: present_oscillator_strengths, present_oscillator_strengths_velocity
      logical :: present_transition_dipoles, present_transition_velocities
      logical :: present_transition_dipole_origin, present_nto_leading_weight
      logical :: present_state_spin
         !! Presence bytes for the excited-state arrays, read off the
         !! wire ahead of each one. A spectrum can arrive without its
         !! transition properties; see the note at the send.

      ! Start from an empty result. pic-mpi's array receives allocate only when
      ! the target is unallocated -- they will not resize a buffer that is
      ! already there, yet still receive the sender's full element count into
      ! it, so reusing one result across messages of different shapes would
      ! overflow the smaller allocation and corrupt the heap.
      call result%destroy()

      ! Receive SCF energy (non-blocking)
      call irecv(comm, result%energy%scf, source, tag, req)

      ! Receive other energy components (blocking to avoid needing multiple request handles)
      call recv(comm, result%energy%dispersion, source, tag, status)
      call recv(comm, result%energy%mp2%ss, source, tag, status)
      call recv(comm, result%energy%mp2%os, source, tag, status)
      call recv(comm, result%energy%cc%singles, source, tag, status)
      call recv(comm, result%energy%cc%doubles, source, tag, status)
      call recv(comm, result%energy%cc%triples, source, tag, status)
      result%has_energy = .true.

      ! Receive fragment metadata
      call recv(comm, result%distance, source, tag, status)
      call recv(comm, result%scf_status, source, tag, status)
      call recv(comm, result%scf_iterations, source, tag, status)
      call recv(comm, result%has_orbitals, source, tag, status)
      call recv(comm, result%homo, source, tag, status)
      call recv(comm, result%lumo, source, tag, status)

      ! Receive gradient flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_gradient, source, tag, status)
      if (result%has_gradient) then
         ! Receive allocatable gradient array (MPI lib handles allocation)
         call recv(comm, result%gradient, source, tag, status)
      end if

      ! Receive Hessian flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_hessian, source, tag, status)
      if (result%has_hessian) then
         ! Receive allocatable Hessian array (MPI lib handles allocation)
         call recv(comm, result%hessian, source, tag, status)
      end if

      ! Receive dipole flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_dipole, source, tag, status)
      if (result%has_dipole) then
         ! Receive allocatable dipole array (MPI lib handles allocation)
         call recv(comm, result%dipole, source, tag, status)
      end if

      ! Receive dipole derivatives flag and data (blocking to avoid needing multiple request handles)
      call recv(comm, result%has_dipole_derivatives, source, tag, status)
      if (result%has_dipole_derivatives) then
         ! Receive allocatable dipole derivatives array (MPI lib handles allocation)
         call recv(comm, result%dipole_derivatives, source, tag, status)
      end if

      ! Receive excited states, every array under one flag
      call recv(comm, result%has_excited_states, source, tag, status)
      if (result%has_excited_states) then
         call recv(comm, present_excitation_energies, source, tag, status)
         if (present_excitation_energies) then
            call recv(comm, result%excitation_energies, source, tag, status)
         end if
         call recv(comm, present_excited_total_energies, source, tag, status)
         if (present_excited_total_energies) then
            call recv(comm, result%excited_total_energies, source, tag, status)
         end if
         call recv(comm, present_oscillator_strengths, source, tag, status)
         if (present_oscillator_strengths) then
            call recv(comm, result%oscillator_strengths, source, tag, status)
         end if
         call recv(comm, present_oscillator_strengths_velocity, source, tag, status)
         if (present_oscillator_strengths_velocity) then
            call recv(comm, result%oscillator_strengths_velocity, source, tag, status)
         end if
         call recv(comm, present_transition_dipoles, source, tag, status)
         if (present_transition_dipoles) then
            call recv(comm, result%transition_dipoles, source, tag, status)
         end if
         call recv(comm, present_transition_velocities, source, tag, status)
         if (present_transition_velocities) then
            call recv(comm, result%transition_velocities, source, tag, status)
         end if
         call recv(comm, present_transition_dipole_origin, source, tag, status)
         if (present_transition_dipole_origin) then
            call recv(comm, result%transition_dipole_origin, source, tag, status)
         end if
         call recv(comm, present_nto_leading_weight, source, tag, status)
         if (present_nto_leading_weight) then
            call recv(comm, result%nto_leading_weight, source, tag, status)
         end if
         call recv(comm, present_state_spin, source, tag, status)
         if (present_state_spin) then
            call recv(comm, result%state_spin, source, tag, status)
         end if
      end if

      call recv_error_state(result, comm, source, tag)
   end subroutine result_irecv

   pure function scf_not_converged_message(iterations) result(text)
      !! What to say when an SCF runs out of cycles and the run must stop
      !!
      !! Shared by every backend. **The escape hatch is named, and so is its
      !! cost**: `allow_crap_scf` lets a fragmented run finish, and the total it
      !! then reports is built partly from unconverged fragments.
      integer, intent(in) :: iterations
      character(len=:), allocatable :: text

      character(len=16) :: cycles

      write (cycles, "(i0)") iterations
      text = "the SCF did not converge in "//trim(cycles)//" cycles at the "// &
             "requested threshold. Raise keywords.scf.maxiter, or loosen "// &
             "keywords.scf.tolerance, or -- in a fragmented run where a few "// &
             "fragments out of many will not converge whatever you do -- set "// &
             "keywords.scf.allow_crap_scf to true. That lets the calculation "// &
             "finish and records every fragment that failed in the output "// &
             "under 'unconverged', with the monomers each was built from and "// &
             "how much energy they carry, which is enough for the Python "// &
             "interface to build a follow-up job over just those. The total "// &
             "reported by such a run is built partly from unconverged "// &
             "fragments, so it is a result to follow up rather than to trust."
   end function scf_not_converged_message

   pure function scf_status_label(status) result(text)
      !! One word for a status, for a table column or a log line
      integer, intent(in) :: status
      character(len=3) :: text

      select case (status)
      case (SCF_CONVERGED)
         text = "yes"
      case (SCF_NOT_CONVERGED)
         text = "NO "
      case default
         text = "?  "
      end select
   end function scf_status_label

   pure subroutine frontier_orbitals(energies, occupations, homo, lumo, found)
      !! Pick HOMO and LUMO out of an ascending orbital spectrum
      !!
      !! Occupation-driven rather than counting electrons, so the same routine
      !! serves restricted and unrestricted spectra. The threshold is half an
      !! electron, which under Fermi smearing names the same pair as any other
      !! cut for a molecule with a real gap, and is arbitrary without one.
      real(dp), intent(in) :: energies(:)
      real(dp), intent(in) :: occupations(:)
      real(dp), intent(out) :: homo, lumo
      logical, intent(out) :: found

      integer :: i, i_homo

      homo = 0.0_dp
      lumo = 0.0_dp
      found = .false.
      i_homo = 0
      do i = 1, min(size(energies), size(occupations))
         if (occupations(i) > 0.5_dp) i_homo = i
      end do

      ! Nothing occupied, or nothing empty: no frontier pair to report.
      if (i_homo < 1 .or. i_homo >= size(energies)) return

      homo = energies(i_homo)
      lumo = energies(i_homo + 1)
      found = .true.
   end subroutine frontier_orbitals

end module mqc_result_types
