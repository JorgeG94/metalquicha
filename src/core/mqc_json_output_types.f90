!! Unified JSON output data container
!! Centralizes all data needed for JSON output from any calculation type
module mqc_json_output_types
   use pic_types, only: int64, dp
   use mqc_thermochemistry, only: thermochemistry_result_t
   implicit none
   private

   public :: json_output_data_t
   public :: OUTPUT_MODE_NONE, OUTPUT_MODE_UNFRAGMENTED, OUTPUT_MODE_MBE, OUTPUT_MODE_GMBE_PIE

   ! Output mode constants
   integer, parameter :: OUTPUT_MODE_NONE = 0
   integer, parameter :: OUTPUT_MODE_UNFRAGMENTED = 1
   integer, parameter :: OUTPUT_MODE_MBE = 2
   integer, parameter :: OUTPUT_MODE_GMBE_PIE = 3

   type :: json_output_data_t
      !! Unified container for all JSON output data
      !!
      !! This type consolidates all data needed to write JSON output for any
      !! calculation type (unfragmented, MBE, GMBE PIE). The output_mode field
      !! determines which format to use when writing JSON.

      integer :: output_mode = OUTPUT_MODE_NONE  !! OUTPUT_MODE_* constant

      !----- Common data -----
      real(dp) :: total_energy = 0.0_dp
      real(dp), allocatable :: gradient(:, :)     !! (3, natoms)
      real(dp), allocatable :: hessian(:, :)      !! (3*natoms, 3*natoms)
      real(dp), allocatable :: dipole(:)          !! (3)
      logical :: has_energy = .false.
      logical :: has_gradient = .false.
      logical :: has_hessian = .false.
      logical :: has_dipole = .false.

      !----- Vibrational data (optional) -----
      real(dp), allocatable :: frequencies(:)        !! cm^-1
      real(dp), allocatable :: reduced_masses(:)     !! amu
      real(dp), allocatable :: force_constants(:)    !! mdyne/Angstrom
      real(dp), allocatable :: ir_intensities(:)     !! km/mol
      type(thermochemistry_result_t) :: thermo
      logical :: has_vibrational = .false.
      logical :: has_ir_intensities = .false.

      !----- MBE-specific data (store ALL fragments for detailed output) -----
      integer, allocatable :: polymers(:, :)          !! Fragment composition (n_fragments, max_level)
      real(dp), allocatable :: fragment_energies(:)   !! Per-fragment total energies
      real(dp), allocatable :: delta_energies(:)      !! MBE delta corrections
      real(dp), allocatable :: sum_by_level(:)        !! Energy sum per level
      real(dp), allocatable :: fragment_distances(:)  !! Per-fragment min distances (Angstrom)
      integer, allocatable :: fragment_charges(:)         !! Per-fragment total charge
      integer, allocatable :: fragment_multiplicities(:)  !! Per-fragment spin multiplicity
      real(dp) :: homo = 0.0_dp          !! Whole-system HOMO, unfragmented runs only
      real(dp) :: lumo = 0.0_dp          !! Whole-system LUMO, unfragmented runs only
      logical :: has_orbitals = .false.
         !! Set only where a gap means something -- one SCF over one system.
         !! A fragmented run leaves this false on purpose: gaps do not add,
         !! so there is no expansion of them to report.
      logical, allocatable :: fragment_has_orbitals(:)
         !! Whether that fragment reported a frontier pair. Not inferred from
         !! the values: homo == lumo == 0 is what a method that said nothing
         !! leaves behind, and printing it as a gap of zero is a claim.
      real(dp), allocatable :: fragment_homo(:)   !! Per-fragment HOMO (Hartree)
      real(dp), allocatable :: fragment_lumo(:)   !! Per-fragment LUMO (Hartree)
      integer, allocatable :: fragment_scf_status(:)
         !! Per-fragment SCF convergence, as SCF_* from mqc_result_types.
         !! Recorded because a non-converged fragment still yields a number of
         !! the right magnitude, so nothing downstream can tell -- and at
         !! millions of terms nobody is reading the log.
      integer(int64) :: fragment_count = 0
      integer :: max_level = 0
      character(len=16) :: fragment_breakdown = "csv"
      character(len=16) :: fingerprint = ""
         !! Identity of the calculation that produced this output. Stamped so a
         !! restart can check what it is about to reuse -- see mqc_fingerprint.
         !! Empty when nothing computed it.
         !! Where the per-fragment table goes: "csv", "json" or "none"

      !----- GMBE PIE-specific data -----
      integer, allocatable :: pie_atom_sets(:, :)     !! Unique atom sets (max_atoms, n_terms)
      integer, allocatable :: pie_coefficients(:)     !! PIE coefficients
      real(dp), allocatable :: pie_energies(:)        !! Per-term energies
      integer(int64) :: n_pie_terms = 0

      !! ---- SAPT ----
      !! An interaction energy decomposed, ordered by `SAPT_TERM_NAMES`. The
      !! total also goes to `total_energy` like any other method's, but on its
      !! own it is the one number a supermolecular calculation would also give;
      !! the decomposition is what the method was run for.
      ! Intrinsic energy decomposition, unfragmented runs that asked for one.
      ! Hartree, and the pair matrices carry the full pair energy in both
      ! (A,B) and (B,A) -- see `calculation_result_t`, which these are copied
      ! from unchanged.
      real(dp), allocatable :: ieda_atom(:)
      real(dp), allocatable :: ieda_free_atom(:)
      real(dp), allocatable :: ieda_pair(:, :)
      real(dp), allocatable :: ieda_classical(:, :)
      real(dp) :: ieda_formation = 0.0_dp
      logical :: has_ieda = .false.

      real(dp), allocatable :: sapt_terms(:)
      logical :: has_sapt = .false.

   contains
      procedure :: destroy => json_output_data_destroy
      procedure :: reset => json_output_data_reset
   end type json_output_data_t

contains

   subroutine json_output_data_destroy(this)
      !! Clean up all allocated memory
      class(json_output_data_t), intent(inout) :: this

      ! Common data
      if (allocated(this%gradient)) deallocate (this%gradient)
      if (allocated(this%hessian)) deallocate (this%hessian)
      if (allocated(this%dipole)) deallocate (this%dipole)

      ! Vibrational data
      if (allocated(this%frequencies)) deallocate (this%frequencies)
      if (allocated(this%reduced_masses)) deallocate (this%reduced_masses)
      if (allocated(this%force_constants)) deallocate (this%force_constants)
      if (allocated(this%ir_intensities)) deallocate (this%ir_intensities)

      ! MBE data
      if (allocated(this%polymers)) deallocate (this%polymers)
      if (allocated(this%fragment_energies)) deallocate (this%fragment_energies)
      if (allocated(this%delta_energies)) deallocate (this%delta_energies)
      if (allocated(this%sum_by_level)) deallocate (this%sum_by_level)
      if (allocated(this%fragment_distances)) deallocate (this%fragment_distances)
      if (allocated(this%fragment_charges)) deallocate (this%fragment_charges)
      if (allocated(this%fragment_multiplicities)) deallocate (this%fragment_multiplicities)

      ! GMBE PIE data
      if (allocated(this%pie_atom_sets)) deallocate (this%pie_atom_sets)
      if (allocated(this%pie_coefficients)) deallocate (this%pie_coefficients)
      if (allocated(this%pie_energies)) deallocate (this%pie_energies)
      if (allocated(this%sapt_terms)) deallocate (this%sapt_terms)
      if (allocated(this%ieda_atom)) deallocate (this%ieda_atom)
      if (allocated(this%ieda_free_atom)) deallocate (this%ieda_free_atom)
      if (allocated(this%ieda_pair)) deallocate (this%ieda_pair)
      if (allocated(this%ieda_classical)) deallocate (this%ieda_classical)

      call this%reset()
   end subroutine json_output_data_destroy

   subroutine json_output_data_reset(this)
      !! Reset all flags and scalar values to defaults
      class(json_output_data_t), intent(inout) :: this

      this%output_mode = OUTPUT_MODE_NONE
      this%total_energy = 0.0_dp
      this%has_energy = .false.
      this%has_gradient = .false.
      this%has_hessian = .false.
      this%has_dipole = .false.
      this%has_vibrational = .false.
      this%has_ir_intensities = .false.
      this%fragment_count = 0
      this%max_level = 0
      this%n_pie_terms = 0
      this%has_sapt = .false.
      this%has_ieda = .false.
      this%ieda_formation = 0.0_dp
   end subroutine json_output_data_reset

end module mqc_json_output_types
