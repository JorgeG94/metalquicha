!! Unified JSON output data container
!! Centralizes all data needed for JSON output from any calculation type
module mqc_json_output_types
   use pic_types, only: int64, dp
   use mqc_thermochemistry, only: thermochemistry_result_t
   use mqc_quao_rows, only: quao_rows_t
   implicit none
   private

   public :: json_output_data_t
   public :: ordered_rows_for_level
   public :: interaction_bonding_term_t
   public :: OUTPUT_MODE_NONE, OUTPUT_MODE_UNFRAGMENTED, OUTPUT_MODE_MBE, OUTPUT_MODE_GMBE_PIE

   ! Output mode constants
   integer, parameter :: OUTPUT_MODE_NONE = 0
   integer, parameter :: OUTPUT_MODE_UNFRAGMENTED = 1
   integer, parameter :: OUTPUT_MODE_MBE = 2
   integer, parameter :: OUTPUT_MODE_GMBE_PIE = 3

   type :: interaction_bonding_term_t
      !! What the bonding analysis of one term says about the reference
      !! fragment's contact with the rest of that term
      !!
      !! Only what crosses between the reference and the other monomers: a row
      !! or pair with both ends on one side is left out, and so is one with an
      !! end on a hydrogen cap or a ghost atom. Atom indices are the system's,
      !! 1-based.
      integer :: term = 0
         !! Row of the term list this came from, 1-based
      integer, allocatable :: monomers(:)
         !! The term's monomers, 1-based
      integer, allocatable :: monomer_of_atom(:, :)
         !! (2, rows%n), the monomer each end of each orbital pair belongs to
      type(quao_rows_t) :: rows
         !! Orbital pairs, `atom` and `partner_atom` renumbered to the system;
         !! `partner_atom` is 0 when the partner is a cap. The atom-pair
         !! matrices are not kept here; see `pair_atoms`.
      integer, allocatable :: pair_atoms(:, :)
         !! (2, n_pairs), reference atom first
      real(dp), allocatable :: pair_bond_index(:)
         !! (n_pairs), sum of squared bond orders between the two atoms
      real(dp), allocatable :: pair_kinetic_bond_order(:)
         !! (n_pairs), kcal/mol
      integer :: omitted_rows = 0
         !! Orbital pairs crossing to the environment that were dropped for
         !! touching a cap or a ghost
   end type interaction_bonding_term_t

   type :: json_output_data_t
      !! Unified container for all JSON output data
      !!
      !! `output_mode` selects the format the writer uses: unfragmented, MBE or
      !! GMBE PIE.

      integer :: output_mode = OUTPUT_MODE_NONE  !! OUTPUT_MODE_* constant

      !----- Common data -----
      real(dp) :: total_energy = 0.0_dp
      real(dp) :: dispersion_energy = 0.0_dp
         !! The empirical dispersion correction, when `keywords.dft.dispersion`
         !! asked for one. Already inside `total_energy`; written separately
         !! because a correction that exists only inside a total cannot be
         !! compared with a published DFT-D number, and a later reader cannot
         !! tell whether a run had one at all.
      logical :: has_dispersion = .false.
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

      !----- Linear-response excited states (optional) -----
      ! Copied unchanged from `calculation_result_t`, which documents the
      ! units, the ordering and the spin codes. Every per-state array runs
      ! over the same states in the same order; the moment arrays are absent
      ! together when the properties could not be computed.
      real(dp), allocatable :: excitation_energies(:)   !! (n_states) Hartree
      real(dp), allocatable :: excited_total_energies(:)  !! (n_states) Hartree
      real(dp), allocatable :: oscillator_strengths(:)  !! (n_states) length gauge
      real(dp), allocatable :: oscillator_strengths_velocity(:)  !! (n_states)
      real(dp), allocatable :: transition_dipoles(:, :)  !! (3, n_states) a.u.
      real(dp), allocatable :: transition_velocities(:, :)  !! (3, n_states) a.u.
      real(dp), allocatable :: transition_dipole_origin(:)  !! (3) Bohr
      real(dp), allocatable :: nto_leading_weight(:)    !! (n_states) 0 to 1
      integer, allocatable :: state_spin(:)             !! (n_states) STATE_SPIN_*
      character(len=16) :: excited_method = ""
         !! Which response problem produced them: "tda" or "rpa". Written
         !! beside the numbers because the two are different answers to the
         !! same deck and nothing in an excitation energy says which.
      character(len=16) :: excited_spin = ""
         !! What the deck asked for: "singlet", "triplet" or "both". The
         !! per-state codes say what each root is; this says what was
         !! requested, which is what a consumer needs to know a list is
         !! complete.
      logical :: has_excited_states = .false.

      !----- MBE-specific data (store ALL fragments for detailed output) -----
      integer, allocatable :: polymers(:, :)          !! Fragment composition (n_fragments, max_level)
      real(dp), allocatable :: fragment_energies(:)   !! Per-fragment total energies
      real(dp), allocatable :: delta_energies(:)      !! MBE delta corrections
      logical, allocatable :: fragment_connected(:)
         !! Whether this term's two monomers are joined by a severed covalent
         !! bond. **Meaningful on two-body rows only** -- false everywhere else,
         !! because the many-body subtraction removes the pair terms and with
         !! them the bond energy. A true row's `delta_energy` includes the
         !! energy of re-forming that bond and is not an interaction energy.
      real(dp), allocatable :: sum_by_level(:)        !! Energy sum per level
      real(dp), allocatable :: fragment_distances(:)  !! Per-fragment min distances (Angstrom)
      integer, allocatable :: fragment_charges(:)         !! Per-fragment total charge
      integer, allocatable :: fragment_multiplicities(:)  !! Per-fragment spin multiplicity
      real(dp) :: homo = 0.0_dp          !! Whole-system HOMO, unfragmented runs only
      real(dp) :: lumo = 0.0_dp          !! Whole-system LUMO, unfragmented runs only
      logical :: has_orbitals = .false.
         !! Set only where a gap means something -- one SCF over one system. A
         !! fragmented run leaves it false: gaps do not add.
      logical, allocatable :: fragment_has_orbitals(:)
         !! Whether that fragment reported a frontier pair. Not inferred from
         !! the values: homo == lumo == 0 is what a method that said nothing
         !! leaves behind, and printing it as a gap of zero is a claim.
      real(dp), allocatable :: fragment_homo(:)   !! Per-fragment HOMO (Hartree)
      real(dp), allocatable :: fragment_lumo(:)   !! Per-fragment LUMO (Hartree)
      integer, allocatable :: fragment_scf_status(:)
         !! Per-fragment SCF convergence, as `SCF_*` from `mqc_result_types`. A
         !! non-converged fragment still yields a number of the right
         !! magnitude, so nothing downstream can tell without this.
      integer(int64), allocatable :: unconverged_ids(:)
         !! Fragment indices whose SCF did not converge, in order, so a
         !! follow-up run can be built from them without reading back a
         !! per-fragment table of millions of rows. Empty when everything
         !! converged, unallocated when the method does not report convergence
         !! at all, which is not the same thing.
      integer, allocatable :: unconverged_monomers(:, :)
         !! (n_unconverged, max_level) the monomers each of those fragments is
         !! built from, zero-padded, exactly as `polymers` holds them.
      real(dp), allocatable :: unconverged_deltas(:)
         !! What each failed fragment contributes to the total, in the same
         !! units and sign as `delta_energies`. The list of failures says which
         !! fragments are suspect; this says whether it matters.
      integer, allocatable :: culprit_monomers(:)
         !! Monomers appearing in at least one failed fragment, most frequent
         !! first, paired with `culprit_counts`.
      integer(int64), allocatable :: culprit_counts(:)
         !! How many failed fragments each of those monomers appears in. **This
         !! collapses a failure list into a diagnosis**: four hundred failures
         !! sharing one monomer is one problem rather than four hundred.
      integer(int64) :: fragment_count = 0
      integer :: max_level = 0

      !----- Interaction energy of one fragment (driver InteractionEnergy) -----
      ! Written in place of `total_energy`, never beside it: the expansion was
      ! reduced to the terms these need, so its sum is not the system's energy
      ! and `has_energy` stays false. See `mbe_result_t`, which these are
      ! copied from.
      logical :: has_interaction = .false.
      integer :: reference_fragment = 0
         !! The reference as a monomer number, 1-based as `polymers` holds it
      real(dp) :: reference_energy = 0.0_dp
      real(dp) :: interaction_energy = 0.0_dp
      real(dp), allocatable :: interaction_by_level(:)       !! (max_level)
      integer(int64), allocatable :: interaction_count_by_level(:)  !! (max_level)
      integer(int64) :: full_expansion_count = 0
         !! How many terms the ordinary expansion would have computed over the
         !! same fragments, level and screening. 0 when not known.
      type(interaction_bonding_term_t), allocatable :: interaction_bonding(:)
         !! One per term holding the reference and at least one other monomer
         !! whose calculation returned a bonding analysis
      logical :: has_interaction_bonding = .false.
      integer, allocatable :: atomic_numbers(:)
         !! (total_atoms), the system's; allocated only when a section names
         !! atoms and wants their elements
      character(len=16) :: fragment_breakdown = "csv"
         !! Where the per-fragment table goes: "csv", "json" or "none"
      character(len=16) :: fingerprint = ""
         !! Identity of the calculation that produced this output. Stamped so a
         !! restart can check what it is about to reuse -- see `mqc_fingerprint`.
         !! Empty when nothing computed it.

      !----- GMBE PIE-specific data -----
      integer, allocatable :: pie_atom_sets(:, :)     !! Unique atom sets (max_atoms, n_terms)
      integer, allocatable :: pie_coefficients(:)     !! PIE coefficients
      real(dp), allocatable :: pie_energies(:)        !! Per-term energies
      integer(int64) :: n_pie_terms = 0

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
      real(dp), allocatable :: atomic_charges(:)
      real(dp), allocatable :: spin_populations(:)
      character(len=16) :: charge_scheme = ""
      logical :: has_charges = .false.
      ! Bond orders over a converged density, when `properties.bond_orders`
      ! asked for them. The scheme travels with the matrix for the reason
      ! given in `calculation_result_t`: two things called bond orders here
      ! are not the same quantity.
      real(dp), allocatable :: bond_orders(:, :)
      real(dp), allocatable :: bond_order_valences(:)
      character(len=16) :: bond_order_scheme = ""
      logical :: has_bond_orders = .false.
      real(dp), allocatable :: fukui_plus(:), fukui_minus(:), fukui_dual(:)
      real(dp) :: fukui_ip = 0.0_dp
      real(dp) :: fukui_ea = 0.0_dp
      real(dp) :: fukui_hardness = 0.0_dp
      real(dp) :: fukui_electrophilicity = 0.0_dp
      logical :: fukui_anion_bound = .true.
      character(len=16) :: fukui_scheme = ""
      logical :: has_fukui = .false.

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

      real(dp), allocatable :: sapt_terms(:)
         !! An interaction energy decomposed, ordered by `SAPT_TERM_NAMES`. The
         !! total also goes to `total_energy` like any other method's, but on
         !! its own it is the one number a supermolecular calculation would also
         !! give; the decomposition is what the method was run for.
      logical :: has_sapt = .false.

      real(dp), allocatable :: efmo_terms(:)
         !! An EFMO energy broken into its sums, ordered by `EFMO_TERM_NAMES`.
         !! The total goes to `total_energy` like any other method's; these are
         !! what says where it came from, and what Phase 3 compares against
         !! GAMESS term by term.
      integer :: efmo_qm_dimers = 0
      integer :: efmo_efp_dimers = 0
         !! How the pairs split at `R_cut`. Together they are every pair, so a
         !! deck can check the cutoff did what was intended without recomputing
         !! the separations.
      integer :: efmo_qm_groups = 0
         !! Near groups of two or more fragments -- the SCFs the near half of
         !! the energy cost. Equal to `efmo_qm_dimers` at level two, larger
         !! above it, and the only place the output says how much a raised
         !! `keywords.fragmentation.level` actually enumerated.
      logical :: has_efmo = .false.

      integer, allocatable :: efmo_pair_fragments(:, :)
         !! (2, n_pairs), the two fragments of each pair, **numbered from one**
         !! as the MBE fragment lists are. Atom indices elsewhere in this file
         !! are 0-based; fragment numbers are not, and these are fragments.
      real(dp), allocatable :: efmo_pair_distance(:)
         !! `R_IJ`, the vdW-scaled closest approach that decided the split.
         !! Unitless, so a reader can see which side of `R_cut` a pair fell.
      logical, allocatable :: efmo_pair_qm(:)
         !! True where a dimer SCF ran, false where four effective-fragment
         !! terms stood in for one.
      real(dp), allocatable :: efmo_pair_energy(:)
         !! What each pair contributed to the total, in Hartree. This is the
         !! interaction map the analysis wants; the aggregate sums beside it
         !! say only what the halves came to.
      integer, allocatable :: efmo_fragment_charges(:)
         !! Each fragment's net charge, indexed from one as the pairs are.
         !! Written beside the pair map because a charged fragment is what
         !! explains a monopole-sized electrostatics term on a pair, and the
         !! pair rows alone cannot say which fragments carry one.
      real(dp), allocatable :: efmo_pair_terms(:, :)
         !! (4, n_pairs): electrostatics, dispersion, exchange repulsion and
         !! charge transfer. Meaningful only where `efmo_pair_qm` is false.

   contains
      procedure :: destroy => json_output_data_destroy
      procedure :: reset => json_output_data_reset
   end type json_output_data_t

contains

   subroutine ordered_rows_for_level(data, frag_level, order)
      !! Row indices of one level's terms, strongest interaction first
      !!
      !! Sorted by the magnitude of the many-body correction, descending, so a
      !! table of thousands of rows is read from the top. Nothing else in this
      !! writer sorts, and nothing in this project did before: the enumeration
      !! order the terms arrive in puts the largest of them nowhere in
      !! particular.
      !!
      !! **Two-body terms joined by a severed covalent bond go last**,
      !! whatever their magnitude, and they are the largest rows in any
      !! fragmented peptide. Their correction is dominated by the energy of
      !! re-forming the bond and is not an interaction energy, so letting them
      !! lead a table sorted by strength would put the one row a reader must
      !! not take at face value at the top of it.
      !!
      !! Level one keeps the order it was enumerated in. A monomer has no
      !! correction to sort by, and its own energy is not a strength.
      use pic_types, only: int_index
      use pic_sorting, only: sort_index
      type(json_output_data_t), intent(in) :: data
      integer, intent(in) :: frag_level
      integer(int64), allocatable, intent(out) :: order(:)

      integer(int64), allocatable :: rows(:)
      integer(int_index), allocatable :: perm(:)
      real(dp), allocatable :: key(:)
      logical, allocatable :: late(:)
      integer(int64) :: i
      integer :: n, k, p

      n = 0
      allocate (rows(data%fragment_count))
      do i = 1_int64, data%fragment_count
         if (count(data%polymers(i, :) > 0) /= frag_level) cycle
         n = n + 1
         rows(n) = i
      end do

      if (frag_level < 2 .or. n < 2 .or. .not. allocated(data%delta_energies)) then
         order = rows(1:n)
         return
      end if

      allocate (key(n), perm(n), late(n))
      do k = 1, n
         key(k) = abs(data%delta_energies(rows(k)))
         late(k) = .false.
         if (frag_level == 2 .and. allocated(data%fragment_connected)) then
            late(k) = data%fragment_connected(rows(k))
         end if
      end do
      call sort_index(key, perm, reverse=.true.)

      ! Two passes over one sorted permutation rather than a compound key: the
      ! grouping stays exact however the magnitudes fall, and each group is
      ! still strongest first within itself.
      allocate (order(n))
      p = 0
      do k = 1, n
         if (late(int(perm(k)))) cycle
         p = p + 1
         order(p) = rows(int(perm(k)))
      end do
      do k = 1, n
         if (.not. late(int(perm(k)))) cycle
         p = p + 1
         order(p) = rows(int(perm(k)))
      end do
   end subroutine ordered_rows_for_level

   subroutine json_output_data_destroy(this)
      !! Clean up all allocated memory
      ! TODO(mqc): every per-fragment array added for SCF status is missing
      ! here -- `fragment_has_orbitals`, `fragment_homo`, `fragment_lumo`,
      ! `fragment_scf_status`, the five `unconverged_*`/`culprit_*` arrays --
      ! so on a reused container they survive with the previous run's contents
      ! and length.
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
      if (allocated(this%fragment_connected)) deallocate (this%fragment_connected)
      if (allocated(this%sum_by_level)) deallocate (this%sum_by_level)
      if (allocated(this%fragment_distances)) deallocate (this%fragment_distances)
      if (allocated(this%fragment_charges)) deallocate (this%fragment_charges)
      if (allocated(this%fragment_multiplicities)) deallocate (this%fragment_multiplicities)
      if (allocated(this%interaction_by_level)) deallocate (this%interaction_by_level)
      if (allocated(this%interaction_count_by_level)) deallocate (this%interaction_count_by_level)
      if (allocated(this%interaction_bonding)) deallocate (this%interaction_bonding)
      if (allocated(this%atomic_numbers)) deallocate (this%atomic_numbers)

      ! GMBE PIE data
      if (allocated(this%pie_atom_sets)) deallocate (this%pie_atom_sets)
      if (allocated(this%pie_coefficients)) deallocate (this%pie_coefficients)
      if (allocated(this%pie_energies)) deallocate (this%pie_energies)
      if (allocated(this%sapt_terms)) deallocate (this%sapt_terms)
      if (allocated(this%efmo_terms)) deallocate (this%efmo_terms)
      if (allocated(this%efmo_pair_fragments)) deallocate (this%efmo_pair_fragments)
      if (allocated(this%efmo_pair_distance)) deallocate (this%efmo_pair_distance)
      if (allocated(this%efmo_pair_qm)) deallocate (this%efmo_pair_qm)
      if (allocated(this%efmo_pair_energy)) deallocate (this%efmo_pair_energy)
      if (allocated(this%efmo_pair_terms)) deallocate (this%efmo_pair_terms)
      if (allocated(this%efmo_fragment_charges)) deallocate (this%efmo_fragment_charges)
      if (allocated(this%ieda_atom)) deallocate (this%ieda_atom)
      if (allocated(this%atomic_charges)) deallocate (this%atomic_charges)
      if (allocated(this%spin_populations)) deallocate (this%spin_populations)
      if (allocated(this%bond_orders)) deallocate (this%bond_orders)
      if (allocated(this%bond_order_valences)) deallocate (this%bond_order_valences)
      if (allocated(this%fukui_plus)) deallocate (this%fukui_plus)
      if (allocated(this%fukui_minus)) deallocate (this%fukui_minus)
      if (allocated(this%fukui_dual)) deallocate (this%fukui_dual)
      if (allocated(this%ieda_free_atom)) deallocate (this%ieda_free_atom)
      if (allocated(this%ieda_pair)) deallocate (this%ieda_pair)
      if (allocated(this%ieda_classical)) deallocate (this%ieda_classical)
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
   end subroutine json_output_data_destroy

   subroutine json_output_data_reset(this)
      !! Reset all flags and scalar values to defaults
      ! TODO(mqc): `has_charges`, `has_orbitals`, `charge_scheme`,
      ! `fukui_scheme` and `fingerprint` are not among them, so a reused
      ! container reports the previous run's charges and gap as its own.
      class(json_output_data_t), intent(inout) :: this

      this%output_mode = OUTPUT_MODE_NONE
      this%total_energy = 0.0_dp
      this%dispersion_energy = 0.0_dp
      this%has_dispersion = .false.
      this%has_energy = .false.
      this%has_gradient = .false.
      this%has_hessian = .false.
      this%has_dipole = .false.
      this%has_vibrational = .false.
      this%has_ir_intensities = .false.
      this%has_excited_states = .false.
      this%excited_method = ""
      this%excited_spin = ""
      this%fragment_count = 0
      this%max_level = 0
      this%has_interaction = .false.
      this%has_interaction_bonding = .false.
      this%reference_fragment = 0
      this%reference_energy = 0.0_dp
      this%interaction_energy = 0.0_dp
      this%full_expansion_count = 0
      this%n_pie_terms = 0
      this%has_sapt = .false.
      this%has_efmo = .false.
      this%efmo_qm_dimers = 0
      this%efmo_efp_dimers = 0
      this%efmo_qm_groups = 0
      this%has_ieda = .false.
      this%has_fukui = .false.
      this%has_bond_orders = .false.
      this%has_stability = .false.
      this%stability_stable = .true.
      this%stability_has_curvature = .false.
      this%stability_curvature = 0.0_dp
      this%stability_rotations = 0
      this%ieda_formation = 0.0_dp
   end subroutine json_output_data_reset

end module mqc_json_output_types
