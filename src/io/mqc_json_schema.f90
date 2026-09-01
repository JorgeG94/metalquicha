!! Schema validation for MQC JSON input decks
module mqc_json_schema
   !! Checks a deck against the input schema before anything reads it.
   !!
   !! Why this is separate from the reader. The reader's `optional_*`
   !! accessors leave a component alone when a key is absent, which is what
   !! keeps the defaults declared in one place -- but it means a key the reader
   !! does not recognise is indistinguishable from a key nobody wrote.
   !! `"maxiter"` misspelled as `"max_iter"` would silently run 300 iterations
   !! instead of the 40 you asked for, and nothing would say so. Rejecting
   !! unknown keys is the only way to tell those two apart, and it cannot be
   !! done by the reader, which by construction only ever looks up names it
   !! already knows.
   !!
   !! So this module walks the document instead of querying it. That needs the
   !! `json_core`/`json_value` interface rather than `json_file`: the latter
   !! reaches a value by path and has no way to ask an object what keys it
   !! actually has.
   !!
   !! What is checked:
   !!
   !!   * every key, at every level, is one the schema defines
   !!   * required keys are present
   !!   * the parallel fragment lists have matching lengths
   !!   * fragment charges sum to the molecular charge
   !!
   !! What is not: types and ranges. The reader reports those where it reads
   !! them, with the value in hand, and duplicating the knowledge here would
   !! give two places to update and two chances to disagree.
   use mqc_error, only: error_t, ERROR_VALIDATION
   use json_module, only: json_file, json_core, json_value
   use mqc_string_utils, only: int_to_text
   implicit none
   private

   public :: ensure_valid_json  !! Reject a deck the schema does not describe

   integer, parameter :: MAX_KEYS = 40
      !! Widest allow-list below; raising it is the only cost of a new key.
      !! optimization_keys is now the widest: the transition-state work added
      !! `target`, and the chain-of-states work `endpoint`, `neb_endpoints`,
      !! `images`, `neb_spring`, the dimer controls and `connect`, which took it
      !! past fragmentation_keys (21,
      !! pushed there by its FMO inner-SCF controls and bond_breaking/cap_scale).
      !! This keeps headroom above that.
      !!
      !! `allow` checks this rather than trusting the comment. It used not to,
      !! and going one over wrote past `allowed` and corrupted whatever followed
      !! -- which showed up as a segfault in an unrelated test four cases later,
      !! not as anything pointing at the key that overflowed.

   type :: key_set_t
      !! The keys one object may contain, and which of them it must
      character(len=32) :: allowed(MAX_KEYS) = ""
      integer :: n_allowed = 0
      character(len=32) :: required(MAX_KEYS) = ""
      integer :: n_required = 0
   end type key_set_t

contains

   subroutine ensure_valid_json(json, error, settings_only)
      !! Validate a loaded deck, setting `error` on the first problem found
      !!
      !! Reports one problem rather than all of them: the first unknown key is
      !! usually the whole story, and a reader that stops there gives a
      !! message short enough to act on.
      !!
      !! `settings_only` validates a document that describes what to compute
      !! but not what to compute it on -- the form that crosses to the workers
      !! when the system arrives by another route. It relaxes exactly one
      !! thing, the requirement that molecules be present; every other key is
      !! held to the same spelling and the same shape, which is the point of
      !! sending a document rather than a struct.
      type(json_file), intent(inout) :: json
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: settings_only

      type(json_core) :: core
      type(json_value), pointer :: root
      logical :: settings

      settings = .false.
      if (present(settings_only)) settings = settings_only

      call json%get_core(core)
      call json%get(root)
      if (.not. associated(root)) then
         call error%set(ERROR_VALIDATION, "Input file contains no JSON document")
         return
      end if

      call check_object(core, root, "", root_keys(settings), error)
      if (error%has_error()) return

      call check_child_object(core, root, "schema", schema_keys(), error)
      if (error%has_error()) return
      call check_child_object(core, root, "model", model_keys(), error)
      if (error%has_error()) return
      call check_child_object(core, root, "system", system_keys(), error)

      call check_child_object(core, root, "properties", properties_keys(), error)
      call check_grandchild_object(core, root, "properties", "fukui", &
                                   fukui_keys(), error)
      if (error%has_error()) return
      call check_fukui_scf_object(core, root, error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "properties", "charges", &
                                   charges_keys(), error)
      call check_grandchild_object(core, root, "properties", "bonding_analysis", &
                                   bonding_analysis_keys(), error)
      call check_bonding_analysis(core, root, error)
      call check_fukui_guess(core, root, error)
      call check_ecp_supported(core, root, error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "system", "logger", logger_keys(), error)
      if (error%has_error()) return

      call check_child_object(core, root, "keywords", keywords_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "scf", scf_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "correlation", &
                                   correlation_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "efp", efp_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "dft", dft_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "pcm", pcm_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "guess", guess_keys(), error)
      if (error%has_error()) return
      call validate_guess_group(core, root, error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "cc", cc_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "mcscf", mcscf_keys(), error)
      if (error%has_error()) return
      call check_avas_group(core, root, error)
      call check_ormas_group(core, root, error)
      call check_valence_choice(core, root, error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "hessian", hessian_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "aimd", aimd_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "optimization", &
                                   optimization_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "xtb", xtb_keys(), error)
      if (error%has_error()) return
      call check_grandchild_object(core, root, "keywords", "fragmentation", &
                                   fragmentation_keys(), error)
      if (error%has_error()) return
      call check_cutoffs(core, root, error)
      if (error%has_error()) return

      call check_molecules(core, root, error)
   end subroutine ensure_valid_json

   ! ==========================================================================
   !  The schema
   !
   !  One function per object, so adding a key is a one-line change next to
   !  the other keys of the same object rather than an edit to a table
   !  somewhere else.
   ! ==========================================================================

   function root_keys(settings_only) result(keys)
      logical, intent(in) :: settings_only
      type(key_set_t) :: keys
      call allow(keys, "schema")
      call allow(keys, "model")
      call allow(keys, "driver")
      call allow(keys, "backend")
      call allow(keys, "molecules")
      call allow(keys, "keywords")
      call allow(keys, "properties")
      call allow(keys, "system")
      call allow(keys, "title")
      call require(keys, "schema")
      call require(keys, "model")
      call require(keys, "driver")
      ! Still *allowed* when settings-only, so a document that carries one gets
      ! the reader's explanation of why it will not be used rather than a bare
      ! "unknown key: molecules", which would read as the opposite problem.
      if (.not. settings_only) call require(keys, "molecules")
   end function root_keys

   function schema_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "name")
      call allow(keys, "version")
      call allow(keys, "index_base")
      call allow(keys, "units")
      call require(keys, "name")
      call require(keys, "version")
   end function schema_keys

   function model_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "method")
      call allow(keys, "basis")
      call allow(keys, "ecp")
      call allow(keys, "aux_basis")
      call allow(keys, "functional")
      call allow(keys, "cartesian")
      call allow(keys, "unrestricted")
         !! Here rather than under `keywords.scf`, where it used to live. That
         !! block is how an answer is reached -- iterations, tolerances, the
         !! extrapolation -- and none of it changes the answer. This does: an
         !! unrestricted reference is a larger variational space and a different
         !! wavefunction, so a broken-symmetry singlet has its own energy. It
         !! belongs with `method` and `functional`, which is what the code has
         !! always said -- `run_libcint_rhf` and `run_libcint_uhf` are separate
         !! routines and the output banner labels runs RKS or UKS.
      call require(keys, "method")
   end function model_keys

   function system_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "logger")
      call allow(keys, "gpu")
      call allow(keys, "skip_json_output")
      call allow(keys, "unchecked_input")
      call allow(keys, "checkpoint")
      call allow(keys, "fragment_breakdown")
   end function system_keys

   function logger_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "level")
   end function logger_keys

   function keywords_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "scf")
      call allow(keys, "hessian")
      call allow(keys, "aimd")
      call allow(keys, "optimization")
      call allow(keys, "fragmentation")
      call allow(keys, "xtb")
      call allow(keys, "correlation")
      call allow(keys, "cc")
      call allow(keys, "efp")
      call allow(keys, "mcscf")
      call allow(keys, "dft")
      call allow(keys, "pcm")
      call allow(keys, "guess")
   end function keywords_keys

   function properties_keys() result(keys)
      !! Things to report once the wave function exists
      !!
      !! Beside `keywords` rather than inside it, and the distinction is worth
      !! keeping: `keywords` say how to compute the wave function and change the
      !! number that comes out, while `properties` ask for something further to
      !! be done with it and change nothing about the energy. A bonding analysis
      !! belongs in the second group, which is why `driver` stays "energy".
      type(key_set_t) :: keys
      call allow(keys, "bonding_analysis")
      call allow(keys, "fukui")
      call allow(keys, "charges")
   end function properties_keys

   function fukui_keys() result(keys)
      !! Settings for the Fukui index analysis
      !!
      !! `population` defaults to CHELPG. It was required rather than
      !! defaulted, on the argument that a choice which changes the numbers
      !! should be written down rather than inherited -- but the two schemes
      !! are not equal candidates. CHELPG fits the electrostatic potential and
      !! is what a condensed Fukui index is normally reported from; Mulliken is
      !! basis-set dependent to the point of changing which site ranks first.
      !! Defaulting to the better one and naming the other is more useful than
      !! refusing to choose, and the scheme is echoed in the report and the
      !! JSON either way, so nothing is inherited silently.
      type(key_set_t) :: keys
      call allow(keys, "population")
      ! Which starting density each ion gets: "neutral" or "independent".
      !
      ! **Deliberately a sibling of `scf`, not a key inside it**, even though
      ! `scf` has a `guess` of its own. The two are different questions and the
      ! spelling collides: this one says whether the ion starts from the
      ! neutral's converged orbitals at all, while `scf.guess` says how a
      ! density is built when it starts from nothing -- `core`, `gwh`, `sad`.
      ! Nesting is what keeps them apart, so moving this one inside `scf` would
      ! put two unrelated meanings on one name.
      call allow(keys, "guess")
      ! The ions are their own SCF problem and are configured as one, with the
      ! same key spellings `keywords.scf` uses so a block can be moved between
      ! them unchanged.
      call allow(keys, "scf")
   end function fukui_keys

   function fukui_scf_keys() result(keys)
      !! `properties.fukui.scf`: how the deltaSCF ions converge
      !!
      !! Every key here is spelled exactly as its `keywords.scf` counterpart,
      !! and each is *seeded* from that counterpart before the deck is read --
      !! so naming none of them leaves the ions converging the way the neutral
      !! did, and naming one changes only that one. There is no sentinel for
      !! "absent" anywhere below, because absent is spelled by the field still
      !! holding what it inherited.
      !!
      !! `inherit_scf` is the exception, and the only key here without a
      !! `keywords.scf` twin: false drops the seed, so unnamed keys fall back
      !! to the type defaults rather than to the neutral's settings. That is
      !! for a neutral converged tightly, or with an accelerator chosen for it,
      !! where inheriting is the wrong start rather than a cautious one.
      !!
      !! No `basis`, `density_fitting` or `unrestricted`: the density
      !! difference is taken over the neutral's own basis functions by
      !! construction, and an ion disagreeing with the neutral about the
      !! Hamiltonian would not give a Fukui index -- it would give a
      !! difference between two different calculations.
      type(key_set_t) :: keys
      call allow(keys, "inherit_scf")
      call allow(keys, "maxiter")
      call allow(keys, "tolerance")
      call allow(keys, "density_tolerance")
      call allow(keys, "gradient_tolerance")
      call allow(keys, "linear_dependence_threshold")
      call allow(keys, "level_shift")
      call allow(keys, "diis")
      call allow(keys, "diis_size")
      call allow(keys, "accelerator")
      call allow(keys, "incremental_fock")
      call allow(keys, "allow_crap_scf")
      call allow(keys, "guess")
   end function fukui_scf_keys

   function charges_keys() result(keys)
      !! Settings for the atomic partial charges
      !!
      !! `scheme` defaults to Mulliken, which is the opposite of the Fukui
      !! default and for a reason. A condensed Fukui index is a difference of
      !! two charges, where Mulliken's basis-set sensitivity does not cancel
      !! and can change which site ranks first; asked for on its own, a charge
      !! is usually wanted as the cheap population number, and Mulliken costs
      !! one trace against an overlap that already exists. CHELPG costs an ESP
      !! evaluation on a few thousand points. Naming the cheap one and
      !! documenting the better one is the honest default here.
      type(key_set_t) :: keys
      call allow(keys, "scheme")
   end function charges_keys

   function bonding_analysis_keys() result(keys)
      !! Settings for the quasi-atomic bonding analysis
      !!
      !! An object rather than a bare string, and that is worth getting right
      !! before anything depends on it: `"bonding_analysis": "gms_quao"` reads
      !! well until the first setting has to go somewhere, and then the choice
      !! is between a second key that has to stay in step with it and a breaking
      !! change to decks already written. `type` follows `keywords.guess.type`,
      !! which is the same shape -- a named choice with settings beside it.
      type(key_set_t) :: keys
      call allow(keys, "type")
      call allow(keys, "energy_threshold")
      call allow(keys, "energy_decomposition")
      call allow(keys, "no_sharing")
      call allow(keys, "no_sharing_ci")
      call allow(keys, "restrict_localization")
      call require(keys, "type")
   end function bonding_analysis_keys

   function avas_keys() result(keys)
      !! Choosing the active space by atomic orbital character
      !!
      !! `orbitals` is a list of labels like `["N 2s", "N 2p"]`. Required, since
      !! an AVAS block that names no orbitals is asking for a selection based on
      !! nothing.
      type(key_set_t) :: keys
      call allow(keys, "orbitals")
      call allow(keys, "threshold")
      call require(keys, "orbitals")
   end function avas_keys

   function guess_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "type")
      call allow(keys, "subscf")
   end function guess_keys

   function subscf_keys() result(keys)
      !! Only meaningful under `type: basis_set_projection`, which
      !! `validate_guess_group` enforces -- a `subscf` block beside any other
      !! guess would be read, validated and then silently ignored.
      type(key_set_t) :: keys
      call allow(keys, "steps")
   end function subscf_keys

   function guess_step_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "basis")
      call allow(keys, "maxiter")
      call allow(keys, "tolerance")
      call require(keys, "basis")
   end function guess_step_keys

   function scf_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "maxiter")
      call allow(keys, "tolerance")
      call allow(keys, "density_tolerance")
      call allow(keys, "gradient_tolerance")
      call allow(keys, "linear_dependence_threshold")
      call allow(keys, "guess")
      call allow(keys, "allow_crap_scf")
      call allow(keys, "density_fitting")
      call allow(keys, "level_shift")
      call allow(keys, "diis")
      call allow(keys, "diis_size")
      call allow(keys, "accelerator")
      call allow(keys, "incremental_fock")
   end function scf_keys

   function efp_keys() result(keys)
      !! MAKEFP settings, deliberately not under "scf"
      !!
      !! The SCF a potential runs already reads `keywords.scf.tolerance` and
      !! `keywords.scf.density_tolerance`. Spelling either of them again here would
      !! give a deck two ways to set one number and this validator no way to object,
      !! since both would be keys it knows -- the trap `correlation_keys` records
      !! for `aux_basis`. Everything here belongs to a stage after the SCF.
      type(key_set_t) :: keys
      call allow(keys, "dynamic_tolerance")
      call allow(keys, "dynamic_maxiter")
      call allow(keys, "allow_crap_response")
      call allow(keys, "response_batch")
      call allow(keys, "response")
      call allow(keys, "vdw_scale")
   end function efp_keys

   function correlation_keys() result(keys)
      !! Post-Hartree-Fock settings, deliberately not under "scf"
      !!
      !! The reference and the correlation treatment are configured apart
      !! because they can disagree: a density-fitted Hartree-Fock followed by a
      !! conventional MP2 is a combination someone will ask for, and one
      !! `density_fitting` shared between them could not express it.
      !!
      !! The auxiliary basis is deliberately *not* here. A basis set belongs in
      !! `model` beside the orbital basis it fits, which is where `model.aux_basis`
      !! is, and having two places to name one meant a deck could name both and
      !! silently prefer one.
      type(key_set_t) :: keys
      call allow(keys, "freeze_core")
      call allow(keys, "n_frozen_core")
      call allow(keys, "density_fitting")
      call allow(keys, "scs")
      call allow(keys, "scs_ss")
      call allow(keys, "scs_os")
   end function correlation_keys

   function dft_keys() result(keys)
      !! The integration grid, which is the only thing DFT adds to an SCF
      !!
      !! The functional is `model.functional`, not a keyword: it names the method
      !! as much as the basis does. These are the quadrature it is integrated on.
      type(key_set_t) :: keys
      call allow(keys, "grid_level")
      call allow(keys, "nlc_grid_level")
      call allow(keys, "radial_points")
      call allow(keys, "angular_points")
      call allow(keys, "screening_tolerance")
      call allow(keys, "block_size")
   end function dft_keys

   function pcm_keys() result(keys)
      !! The polarizable continuum: the cavity, the solvent, and the charge solve
      !!
      !! Separate from `xtb`, which configures tblite's own CPCM. Two continuum
      !! implementations with different cavities should not share one keyword.
      type(key_set_t) :: keys
      call allow(keys, "method")
      call allow(keys, "dielectric")
      call allow(keys, "angular_points")
      call allow(keys, "radii_scale")
      call allow(keys, "zeta")
      call allow(keys, "tolerance")
      call allow(keys, "max_iter")
   end function pcm_keys

   function cc_keys() result(keys)
      !! Coupled-cluster settings, separate from "correlation"
      !!
      !! `correlation` holds what every post-Hartree-Fock method shares -- the
      !! frozen core, the fitting, the auxiliary basis. These are the ones only an
      !! iterative method has: how long to iterate and how hard, and whether the
      !! triples correction runs. `triples` is here as an override; ordinarily the
      !! method name settles it, since "ccsd" and "ccsd(t)" are separate method
      !! types rather than one type with a flag.
      type(key_set_t) :: keys
      call allow(keys, "maxiter")
      call allow(keys, "tolerance")
      call allow(keys, "triples")
      call allow(keys, "diis")
      call allow(keys, "diis_size")
      call allow(keys, "spin_adapted")
   end function cc_keys

   function ormas_keys() result(keys)
      !! Cutting the active space into subspaces with occupation windows
      !!
      !! All three are required and all three are lists of the same length:
      !! where each subspace starts, and the fewest and most electrons it may
      !! hold. A partition given by halves is not a partition.
      type(key_set_t) :: keys
      call allow(keys, "subspaces")
      call allow(keys, "min_electrons")
      call allow(keys, "max_electrons")
      call require(keys, "subspaces")
      call require(keys, "min_electrons")
      call require(keys, "max_electrons")
   end function ormas_keys

   function mcscf_keys() result(keys)
      !! The active space, and how hard to work at it
      !!
      !! **Only what the backend actually acts on is listed.** `mcscf_config_t`
      !! carries fields for state averaging and for a CASPT2/NEVPT2 correction,
      !! and none of that is implemented -- `run_libcint_casscf` optimises one
      !! state and there is no perturbative step at all. Allowing the keys would
      !! mean a deck could ask for a three-state average, get a ground-state
      !! energy, and find nothing in the output to say so. Left out, the
      !! validator refuses the key and lists what may be written instead, which
      !! is the difference between "not yet" and "quietly ignored".
      !!
      !! `max_micro_iter` and a CI threshold are absent for a duller reason:
      !! neither routine underneath takes them. The macro loop pins its CASCI at
      !! 1e-11 so the orbital gradient it differentiates is not contaminated by
      !! a loose CI, and that is not a knob worth exposing.
      type(key_set_t) :: keys
      call allow(keys, "avas")
      call allow(keys, "ormas")
      call allow(keys, "full_valence")
      call allow(keys, "n_active_electrons")
      call allow(keys, "n_active_orbitals")
      call allow(keys, "n_inactive_orbitals")
      call allow(keys, "optimize_orbitals")
      call allow(keys, "max_macro_iter")
      call allow(keys, "orbital_convergence")
   end function mcscf_keys

   function hessian_keys() result(keys)
      type(key_set_t) :: keys
      ! Both spellings of the displacement are accepted, as the retired JSON
      ! generator accepted them.
      call allow(keys, "finite_difference_displacement")
      call allow(keys, "displacement")
      call allow(keys, "temperature")
      call allow(keys, "pressure")
   end function hessian_keys

   function aimd_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "dt")
      call allow(keys, "timestep")
      call allow(keys, "nsteps")
      call allow(keys, "steps")
      call allow(keys, "initial_temperature")
      call allow(keys, "temperature")
      call allow(keys, "output_frequency")
      call allow(keys, "output_freq")
   end function aimd_keys

   function optimization_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "max_steps")
      call allow(keys, "steps")
      call allow(keys, "gradient_tolerance")
      call allow(keys, "tolerance")
      call allow(keys, "energy_tolerance")
      call allow(keys, "max_step")
      call allow(keys, "coordinates")
      call allow(keys, "coordinate_system")
      call allow(keys, "algorithm")
      call allow(keys, "optimizer")
      call allow(keys, "lbfgs_memory")
      call allow(keys, "print_level")
      call allow(keys, "trajectory")
      call allow(keys, "freeze_terms")
      call allow(keys, "hess_end")
      call allow(keys, "hessian_update")
      call allow(keys, "target")
      call allow(keys, "endpoint")
      call allow(keys, "neb_endpoints")
      call allow(keys, "images")
      call allow(keys, "neb_spring")
      call allow(keys, "saddle_method")
      call allow(keys, "dimer_separation")
      call allow(keys, "dimer_max_rotations")
      call allow(keys, "dimer_rotation_tolerance")
      call allow(keys, "zero_modes")
      call allow(keys, "soft_mode_threshold")
      call allow(keys, "connect")
      call allow(keys, "connect_distort")
      call allow(keys, "frozen_atoms")
      call allow(keys, "constraints")
      call allow(keys, "timestep")
      call allow(keys, "friction")
      call allow(keys, "friction_factor")
      call allow(keys, "friction_rising")
   end function optimization_keys

   function xtb_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "solvent")
      call allow(keys, "solvation_model")
      call allow(keys, "dielectric")
      call allow(keys, "cpcm_nang")
      call allow(keys, "cpcm_rscale")
      call allow(keys, "use_cds")
      call allow(keys, "use_shift")
   end function xtb_keys

   function fragmentation_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "method")
      call allow(keys, "level")
      call allow(keys, "allow_overlapping_fragments")
      call allow(keys, "max_intersection_level")
      call allow(keys, "expansion")
      call allow(keys, "counterpoise")
      call allow(keys, "far_field")
      call allow(keys, "resppc")
      call allow(keys, "max_outer")
      call allow(keys, "outer_tolerance")
      call allow(keys, "scf_max_iter")
      call allow(keys, "scf_energy_tolerance")
      call allow(keys, "scf_density_tolerance")
      call allow(keys, "bond_breaking")
      call allow(keys, "cap_scale")
      call allow(keys, "embedding")
      call allow(keys, "cutoff_method")
      call allow(keys, "distance_metric")
      call allow(keys, "cutoffs")
      call allow(keys, "global_groups")
      call allow(keys, "nodes_per_group")
      call require(keys, "method")
   end function fragmentation_keys

   function molecule_keys() result(keys)
      type(key_set_t) :: keys
      call allow(keys, "symbols")
      call allow(keys, "geometry")
      call allow(keys, "xyz")
      call allow(keys, "molecular_charge")
      call allow(keys, "molecular_multiplicity")
      call allow(keys, "connectivity")
      call allow(keys, "fragments")
      call allow(keys, "fragment_charges")
      call allow(keys, "fragment_multiplicities")
      ! One effective fragment potential per fragment, or a single one standing for all
      ! of them when `uniform_system` says every fragment is the same species -- which
      ! is what keeps a ten-thousand-water deck from repeating one filename ten
      ! thousand times. A fragment with no potential named for it stays quantum, so
      ! this is also how a mixed calculation is spelled.
      call allow(keys, "fragment_potentials")
      call allow(keys, "uniform_system")
      call require(keys, "molecular_charge")
      call require(keys, "molecular_multiplicity")
   end function molecule_keys

   ! ==========================================================================
   !  Walking
   ! ==========================================================================

   subroutine check_object(core, object, path, keys, error)
      !! Every key of `object` must be in `keys`, and every required one present
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: object
      character(len=*), intent(in) :: path  !! For the message; "" at the root
      type(key_set_t), intent(in) :: keys
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: child
      character(len=:), allocatable :: name
      integer :: n_children, ichild, ikey
      logical :: found

      call core%info(object, n_children=n_children)

      do ichild = 1, n_children
         call core%get_child(object, ichild, child, found)
         if (.not. found) cycle
         call core%info(child, name=name)
         if (.not. allocated(name)) cycle
         if (.not. is_allowed(keys, name)) then
            ! A key that moved gets told where it went. The generic message
            ! lists what is allowed *here*, which for a relocated key is the
            ! one list guaranteed not to contain it -- so it reads as "no such
            ! setting" when the setting exists a line away.
            if (path == "keywords.scf" .and. name == "unrestricted") then
               call error%set(ERROR_VALIDATION, "keywords.scf.unrestricted has moved "// &
                              "to model.unrestricted. It selects which wavefunction is "// &
                              "computed, not how the SCF reaches it, so it belongs with "// &
                              "model.method and model.functional.")
               return
            end if
            call error%set(ERROR_VALIDATION, "Unknown key "//quoted(path, name)// &
                           ". Allowed here: "//allowed_list(keys))
            return
         end if
      end do

      do ikey = 1, keys%n_required
         call core%get(object, trim(keys%required(ikey)), child, found)
         if (.not. found .or. .not. associated(child)) then
            call error%set(ERROR_VALIDATION, "Missing required key "// &
                           quoted(path, trim(keys%required(ikey))))
            return
         end if
      end do
   end subroutine check_object

   subroutine check_child_object(core, parent, name, keys, error)
      !! Validate `parent.name` if it is there; absent is fine
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: parent
      character(len=*), intent(in) :: name
      type(key_set_t), intent(in) :: keys
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: child
      logical :: found

      call core%get(parent, name, child, found)
      if (.not. found .or. .not. associated(child)) return
      call check_object(core, child, name, keys, error)
   end subroutine check_child_object

   subroutine check_grandchild_object(core, root, parent_name, name, keys, error)
      !! Validate `parent_name.name` if both levels are there
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      character(len=*), intent(in) :: parent_name, name
      type(key_set_t), intent(in) :: keys
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: parent, child
      logical :: found

      call core%get(root, parent_name, parent, found)
      if (.not. found .or. .not. associated(parent)) return
      call core%get(parent, name, child, found)
      if (.not. found .or. .not. associated(child)) return
      call check_object(core, child, parent_name//"."//name, keys, error)
   end subroutine check_grandchild_object

   subroutine check_avas_group(core, root, error)
      !! `keywords.mcscf.avas`, and that it does not fight the explicit space
      !!
      !! An active space can be named directly, by counts, or chosen by AVAS
      !! from orbital labels. Giving both is refused rather than resolved by
      !! precedence: two ways of saying one thing leaves the deck's meaning
      !! depending on which the reader reaches first, and the counts a user
      !! wrote would be silently discarded in favour of whatever the projection
      !! decided.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: keywords, mcscf, avas, entry
      logical :: found, has_counts

      if (error%has_error()) return
      call core%get(root, "keywords", keywords, found)
      if (.not. found .or. .not. associated(keywords)) return
      call core%get(keywords, "mcscf", mcscf, found)
      if (.not. found .or. .not. associated(mcscf)) return
      call core%get(mcscf, "avas", avas, found)
      if (.not. found .or. .not. associated(avas)) return

      call check_object(core, avas, "keywords.mcscf.avas", avas_keys(), error)
      if (error%has_error()) return

      has_counts = .false.
      call core%get(mcscf, "n_active_electrons", entry, found)
      if (found .and. associated(entry)) has_counts = .true.
      call core%get(mcscf, "n_active_orbitals", entry, found)
      if (found .and. associated(entry)) has_counts = .true.
      if (has_counts) then
         call error%set(ERROR_VALIDATION, "keywords.mcscf gives both an AVAS block "// &
                        "and an explicit active space. AVAS chooses the space from "// &
                        "the orbital labels, so the counts would be discarded. Give "// &
                        "one or the other.")
         return
      end if

      call core%get(mcscf, "full_valence", entry, found)
      if (found .and. associated(entry)) then
         call error%set(ERROR_VALIDATION, "keywords.mcscf gives both an AVAS block "// &
                        "and full_valence. Each picks the whole active space, so one "// &
                        "of them would be thrown away. Give one or the other.")
         return
      end if

      if (child_count(core, avas, "orbitals") <= 0) then
         call error%set(ERROR_VALIDATION, "keywords.mcscf.avas.orbitals is empty. It "// &
                        "lists the atomic orbitals the active space should be built "// &
                        'from, like ["N 2s", "N 2p"].')
      end if
   end subroutine check_avas_group

   subroutine check_valence_choice(core, root, error)
      !! `keywords.mcscf.full_valence`, and that it is the only space named
      !!
      !! The valence shell is a complete description of an active space, so
      !! counts written beside it would be silently discarded. Refused rather
      !! than resolved, for the reason the AVAS check gives.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: keywords, mcscf, entry
      logical :: found, has_counts

      if (error%has_error()) return
      call core%get(root, "keywords", keywords, found)
      if (.not. found .or. .not. associated(keywords)) return
      call core%get(keywords, "mcscf", mcscf, found)
      if (.not. found .or. .not. associated(mcscf)) return
      call core%get(mcscf, "full_valence", entry, found)
      if (.not. found .or. .not. associated(entry)) return

      has_counts = .false.
      call core%get(mcscf, "n_active_electrons", entry, found)
      if (found .and. associated(entry)) has_counts = .true.
      call core%get(mcscf, "n_active_orbitals", entry, found)
      if (found .and. associated(entry)) has_counts = .true.
      if (has_counts) then
         call error%set(ERROR_VALIDATION, "keywords.mcscf gives both full_valence "// &
                        "and an explicit active space. The valence shell already "// &
                        "says how many orbitals and electrons there are, so the "// &
                        "counts would be discarded. Give one or the other.")
      end if
   end subroutine check_valence_choice

   subroutine check_ormas_group(core, root, error)
      !! `keywords.mcscf.ormas`, and that it describes a partition
      !!
      !! The lengths are checked here rather than left to the builder because
      !! the deck is where the mistake is: three lists of different lengths is
      !! a typo, and saying so with the key names beats an error about class
      !! enumeration from four layers down.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: keywords, mcscf, ormas
      logical :: found
      integer :: n_spaces, n_min, n_max

      if (error%has_error()) return
      call core%get(root, "keywords", keywords, found)
      if (.not. found .or. .not. associated(keywords)) return
      call core%get(keywords, "mcscf", mcscf, found)
      if (.not. found .or. .not. associated(mcscf)) return
      call core%get(mcscf, "ormas", ormas, found)
      if (.not. found .or. .not. associated(ormas)) return

      call check_object(core, ormas, "keywords.mcscf.ormas", ormas_keys(), error)
      if (error%has_error()) return

      n_spaces = child_count(core, ormas, "subspaces")
      n_min = child_count(core, ormas, "min_electrons")
      n_max = child_count(core, ormas, "max_electrons")

      if (n_spaces <= 0) then
         call error%set(ERROR_VALIDATION, "keywords.mcscf.ormas.subspaces is empty. "// &
                        "It lists the active orbital each subspace starts at, so the "// &
                        "first entry is 1.")
         return
      end if
      if (n_min /= n_spaces .or. n_max /= n_spaces) then
         call error%set(ERROR_VALIDATION, "keywords.mcscf.ormas names "// &
                        int_to_text(n_spaces)//" subspaces but gives "// &
                        int_to_text(n_min)//" minima and "//int_to_text(n_max)// &
                        " maxima. All three lists describe the same subspaces and "// &
                        "have to be the same length.")
         return
      end if
   end subroutine check_ormas_group

   pure function lowered(text) result(out)
      !! ASCII lowercase, so a method name matches however a deck spells it
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, c

      out = text
      do i = 1, len(text)
         c = iachar(text(i:i))
         if (c >= iachar("A") .and. c <= iachar("Z")) then
            out(i:i) = achar(c + 32)
         end if
      end do
   end function lowered

   subroutine check_fukui_scf_object(core, root, error)
      !! Validate `properties.fukui.scf`, one level deeper than the helpers go
      !!
      !! `check_grandchild_object` walks two levels and this key is three down,
      !! so the walk is spelled out. Worth validating rather than waving
      !! through: an unrecognised key here is a convergence setting the ions
      !! silently did not get, and the run still prints Fukui indices.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: properties, fukui, scf
      logical :: found

      call core%get(root, "properties", properties, found)
      if (.not. found .or. .not. associated(properties)) return
      call core%get(properties, "fukui", fukui, found)
      if (.not. found .or. .not. associated(fukui)) return
      call core%get(fukui, "scf", scf, found)
      if (.not. found .or. .not. associated(scf)) return
      call check_object(core, scf, "properties.fukui.scf", fukui_scf_keys(), error)
   end subroutine check_fukui_scf_object

   subroutine check_fukui_guess(core, root, error)
      !! `properties.fukui.guess` names a mode this code has
      !!
      !! Checked here rather than left to the reader, for the reason the whole
      !! validator exists: the reader compares against "independent" and treats
      !! everything else as the default, so a misspelling would silently seed
      !! from the neutral -- which is the setting the deck was trying to turn
      !! off. A wrong answer produced by a typo in the thing meant to prevent
      !! it is the worst shape this failure could take.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: properties, fukui, entry
      character(len=:), allocatable :: mode
      logical :: found

      if (error%has_error()) return
      call core%get(root, "properties", properties, found)
      if (.not. found .or. .not. associated(properties)) return
      call core%get(properties, "fukui", fukui, found)
      if (.not. found .or. .not. associated(fukui)) return
      call core%get(fukui, "guess", entry, found)
      if (.not. found .or. .not. associated(entry)) return
      call core%get(fukui, "guess", mode)

      select case (lowered(trim(adjustl(mode))))
      case ("neutral", "independent")
         return
      case default
         call error%set(ERROR_VALIDATION, "properties.fukui.guess is '"//trim(mode)// &
                        "'. It is 'neutral' -- seed each ion from the converged "// &
                        "neutral's orbitals, the default -- or 'independent', which "// &
                        "gives each ion the ordinary SCF guess instead.")
      end select
   end subroutine check_fukui_guess

   subroutine check_ecp_supported(core, root, error)
      !! `model.ecp` only reaches the methods that were wired for it
      !!
      !! Checked here rather than left to each method because the failure is
      !! silence: `ecp_set` is threaded into the Hartree-Fock and DFT builders
      !! and nowhere else, so a CASSCF or xTB deck naming a potential today
      !! runs all-electron and reports a number with nothing to say it ignored
      !! half the request. That is the failure mode the ECP work exists to
      !! prevent, so it must not be the one the deck layer has.
      !!
      !! Correlated methods built on Hartree-Fock -- MP2, RI-MP2, coupled
      !! cluster, the double hybrids -- are not listed: they take their
      !! integrals from the SCF, which does carry the potential. Their nuclear
      !! derivatives refuse separately, in the backend.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: model, entry
      character(len=:), allocatable :: ecp, method
      logical :: found

      if (error%has_error()) return
      call core%get(root, "model", model, found)
      if (.not. found .or. .not. associated(model)) return
      call core%get(model, "ecp", entry, found)
      if (.not. found .or. .not. associated(entry)) return
      call core%get(model, "ecp", ecp)
      if (len_trim(ecp) == 0) return

      call core%get(model, "method", entry, found)
      if (.not. found .or. .not. associated(entry)) return
      call core%get(model, "method", method)

      select case (lowered(trim(adjustl(method))))
      case ("gfn1", "gfn2", "gfn0", "xtb")
         call error%set(ERROR_VALIDATION, "model.ecp is set, but '"//trim(method)// &
                        "' is a semiempirical method with its own parameterised "// &
                        "core -- there is no all-electron basis for a potential to "// &
                        "replace. Remove model.ecp, or choose an ab initio method.")
      case ("casscf", "casci", "ormas", "ormas-scf", "ormasscf")
         call error%set(ERROR_VALIDATION, "model.ecp is set, but the MCSCF path does "// &
                        "not carry an effective core potential yet: it would be "// &
                        "ignored and the calculation would run all-electron without "// &
                        "saying so. Refused rather than silently dropped. Hartree-Fock "// &
                        "and DFT, and the correlated methods built on them, do carry it.")
      case default
         ! Everything else carries the potential, or is refused nearer the
         ! calculation: the GPU backend refuses it in its own driver, since
         ! which backend runs is not a `model` key and cannot be seen here.
      end select
   end subroutine check_ecp_supported

   subroutine check_bonding_analysis(core, root, error)
      !! The value of `properties.bonding_analysis` names an analysis we have
      !!
      !! Checked here rather than left to the reader because an unknown name is
      !! the failure this schema exists to catch: a deck that asks for an
      !! analysis nobody implements would otherwise run a perfectly good energy
      !! and report nothing, which looks like the analysis found nothing to say.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: properties, analysis, entry
      character(len=:), allocatable :: name
      logical :: found

      if (error%has_error()) return
      call core%get(root, "properties", properties, found)
      if (.not. found .or. .not. associated(properties)) return
      call core%get(properties, "bonding_analysis", analysis, found)
      if (.not. found .or. .not. associated(analysis)) return
      call core%get(analysis, "type", entry, found)
      if (.not. found .or. .not. associated(entry)) return

      call core%get(analysis, "type", name)
      select case (trim(adjustl(name)))
      case ("none", "gms_quao", "quao")
      case default
         call error%set(ERROR_VALIDATION, "properties.bonding_analysis.type is '"// &
                        trim(name)//"'. Known analyses: none, gms_quao (the "// &
                        "Ruedenberg quasi-atomic bonding picture, spelled as "// &
                        "GAMESS implements it; 'quao' is accepted for it too).")
         return
      end select

      call core%get(analysis, "no_sharing_ci", entry, found)
      if (.not. found .or. .not. associated(entry)) return
      call core%get(analysis, "no_sharing_ci", name)
      select case (trim(adjustl(name)))
      case ("transform", "resolve")
      case default
         call error%set(ERROR_VALIDATION, "properties.bonding_analysis.no_sharing_ci "// &
                        "is '"//trim(name)//"'. It chooses how the CI expansion over "// &
                        "quasi-atomic orbitals is obtained: 'transform' solves in the "// &
                        "molecular orbital basis and carries the vector across with "// &
                        "the orbital transformation, 'resolve' runs a second Davidson "// &
                        "in the quasi-atomic basis. Both describe the same wave "// &
                        "function.")
      end select
   end subroutine check_bonding_analysis

   subroutine validate_guess_group(core, root, error)
      !! The three rules that make `keywords.guess` unambiguous
      !!
      !!   1. `keywords.scf.guess` and `keywords.guess.type` must not both be
      !!      set. The second supersedes the first; refusing both rather than
      !!      picking one keeps the deck's meaning independent of reader order.
      !!   2. `subscf` is only meaningful under `basis_set_projection`. Beside
      !!      any other guess it would be read, validated and then ignored, which
      !!      is the failure mode this whole schema exists to prevent.
      !!   3. `basis_set_projection` needs at least one step. A ladder with no
      !!      rungs is the default guess wearing a longer name.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: keywords, guess, scf, subscf, steps, step
      character(len=:), allocatable :: gtype
      logical :: found, has_scf_guess, is_projection
      integer :: n_steps, i

      call core%get(root, "keywords", keywords, found)
      if (.not. found .or. .not. associated(keywords)) return
      call core%get(keywords, "guess", guess, found)
      if (.not. found .or. .not. associated(guess)) return

      has_scf_guess = .false.
      call core%get(keywords, "scf", scf, found)
      if (found .and. associated(scf)) then
         call core%get(scf, "guess", steps, has_scf_guess)
      end if

      call core%get(guess, "type", steps, found)
      if (found .and. associated(steps)) then
         call core%get(guess, "type", gtype)
      end if

      if (has_scf_guess .and. allocated(gtype)) then
         call error%set(ERROR_VALIDATION, "keywords.scf.guess and keywords.guess.type "// &
                        "both set the initial guess. Use keywords.guess.type; "// &
                        "keywords.scf.guess is the older spelling and is superseded.")
         return
      end if

      is_projection = .false.
      if (allocated(gtype)) is_projection = trim(adjustl(gtype)) == "basis_set_projection"

      call core%get(guess, "subscf", subscf, found)
      if (found .and. associated(subscf)) then
         if (.not. is_projection) then
            call error%set(ERROR_VALIDATION, "keywords.guess.subscf only means something "// &
                           "for type 'basis_set_projection'; beside any other guess it "// &
                           "would be read and then ignored")
            return
         end if
         call check_object(core, subscf, "keywords.guess.subscf", subscf_keys(), error)
         if (error%has_error()) return
         call core%get(subscf, "steps", steps, found)
         if (found .and. associated(steps)) then
            call core%info(steps, n_children=n_steps)
            do i = 1, n_steps
               call core%get_child(steps, i, step, found)
               if (.not. found .or. .not. associated(step)) cycle
               call check_object(core, step, "keywords.guess.subscf.steps", &
                                 guess_step_keys(), error)
               if (error%has_error()) return
            end do
         end if
      end if

      if (is_projection) then
         n_steps = 0
         call core%get(guess, "subscf.steps", steps, found)
         if (found .and. associated(steps)) call core%info(steps, n_children=n_steps)
         if (n_steps < 1) then
            call error%set(ERROR_VALIDATION, "guess type 'basis_set_projection' needs "// &
                           "keywords.guess.subscf.steps with at least one basis to "// &
                           "converge before the target one")
            return
         end if
      end if
   end subroutine validate_guess_group

   subroutine check_cutoffs(core, root, error)
      !! Cutoff keys are n-mer names or decimal levels, so they get their own rule
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: keywords, fragmentation, cutoffs, child
      character(len=:), allocatable :: name
      integer :: n_children, ichild
      logical :: found

      call core%get(root, "keywords", keywords, found)
      if (.not. found .or. .not. associated(keywords)) return
      call core%get(keywords, "fragmentation", fragmentation, found)
      if (.not. found .or. .not. associated(fragmentation)) return
      call core%get(fragmentation, "cutoffs", cutoffs, found)
      if (.not. found .or. .not. associated(cutoffs)) return

      call core%info(cutoffs, n_children=n_children)
      do ichild = 1, n_children
         call core%get_child(cutoffs, ichild, child, found)
         if (.not. found) cycle
         call core%info(child, name=name)
         if (.not. allocated(name)) cycle
         if (nmer_level(name) <= 0) then
            call error%set(ERROR_VALIDATION, "Unknown key "// &
                           quoted("keywords.fragmentation.cutoffs", name)// &
                           ". Use an n-mer name (dimer, trimer, tetramer, "// &
                           "pentamer, hexamer, heptamer, octamer) or a level "// &
                           "number of 2 or more")
            return
         end if
      end do
   end subroutine check_cutoffs

   subroutine check_molecules(core, root, error)
      !! Every molecule's keys, plus the cross-checks between its lists
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: root
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: molecules, molecule
      integer :: n_molecules, imol
      logical :: found

      call core%get(root, "molecules", molecules, found)
      if (.not. found .or. .not. associated(molecules)) return

      call core%info(molecules, n_children=n_molecules)
      if (n_molecules <= 0) then
         call error%set(ERROR_VALIDATION, "'molecules' must contain at least one entry")
         return
      end if

      do imol = 1, n_molecules
         call core%get_child(molecules, imol, molecule, found)
         if (.not. found) cycle
         call check_object(core, molecule, molecule_path(imol), molecule_keys(), error)
         if (error%has_error()) return
         call check_molecule_geometry(core, molecule, imol, error)
         if (error%has_error()) return
         call check_molecule_fragments(core, molecule, imol, error)
         if (error%has_error()) return
      end do
   end subroutine check_molecules

   subroutine check_molecule_geometry(core, molecule, imol, error)
      !! Exactly one of `xyz` or `symbols` + `geometry`
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: molecule
      integer, intent(in) :: imol
      type(error_t), intent(inout) :: error

      logical :: has_xyz, has_symbols, has_geometry

      has_xyz = has_key(core, molecule, "xyz")
      has_symbols = has_key(core, molecule, "symbols")
      has_geometry = has_key(core, molecule, "geometry")

      if (has_xyz .and. (has_symbols .or. has_geometry)) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": give either 'xyz' or 'symbols' with 'geometry', not both")
         return
      end if
      if (.not. has_xyz) then
         if (.not. (has_symbols .and. has_geometry)) then
            call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                           ": must give either 'xyz' or 'symbols' with 'geometry'")
            return
         end if
      end if
   end subroutine check_molecule_geometry

   subroutine check_molecule_fragments(core, molecule, imol, error)
      !! Parallel list lengths, and charges that add up
      !!
      !! Fragment charges summing to something other than the molecular charge
      !! is the one of these that produces a plausible wrong answer rather than
      !! a crash, which is why it is worth catching here.
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: molecule
      integer, intent(in) :: imol
      type(error_t), intent(inout) :: error

      type(json_value), pointer :: fragments, charges
      integer :: n_fragments, n_charges, n_mults, total, molecular_charge
      logical :: found

      call core%get(molecule, "fragments", fragments, found)
      if (.not. found .or. .not. associated(fragments)) return
      call core%info(fragments, n_children=n_fragments)
      if (n_fragments <= 0) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": 'fragments' must not be empty if it is given")
         return
      end if

      n_charges = child_count(core, molecule, "fragment_charges")
      if (n_charges > 0 .and. n_charges /= n_fragments) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": 'fragment_charges' has "//int_to_text(n_charges)// &
                        " entries for "//int_to_text(n_fragments)//" fragments")
         return
      end if

      n_mults = child_count(core, molecule, "fragment_multiplicities")
      if (n_mults > 0 .and. n_mults /= n_fragments) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": 'fragment_multiplicities' has "//int_to_text(n_mults)// &
                        " entries for "//int_to_text(n_fragments)//" fragments")
         return
      end if

      if (n_charges <= 0) return
      call core%get(molecule, "fragment_charges", charges, found)
      if (.not. found) return
      total = sum_integer_array(core, charges, n_charges)

      call core%get(molecule, "molecular_charge", molecular_charge, found)
      if (.not. found) return
      if (total /= molecular_charge) then
         call error%set(ERROR_VALIDATION, molecule_path(imol)// &
                        ": fragment charges sum to "//int_to_text(total)// &
                        " but the molecular charge is "//int_to_text(molecular_charge)// &
                        "; these must agree")
      end if
   end subroutine check_molecule_fragments

   ! ==========================================================================
   !  Helpers
   ! ==========================================================================

   pure subroutine allow(keys, name)
      !! Record a key an object may contain
      !!
      !! Stops rather than overflowing. These sets are built from literals at
      !! startup, so going over MAX_KEYS is a programming error fixed by raising
      !! it -- never something a deck can provoke -- and failing here names the
      !! key that did it instead of corrupting memory for a later test to trip
      !! over.
      type(key_set_t), intent(inout) :: keys
      character(len=*), intent(in) :: name

      if (keys%n_allowed >= MAX_KEYS) then
         error stop "mqc_json_schema: allow-list is full; raise MAX_KEYS"
      end if
      keys%n_allowed = keys%n_allowed + 1
      keys%allowed(keys%n_allowed) = name
   end subroutine allow

   pure subroutine require(keys, name)
      !! Record a key an object must contain
      !!
      !! Required keys are also allowed ones, so every call here is paired with
      !! an `allow` above rather than implying it -- the allow-list is what the
      !! error message offers the user, and it should read as the whole set.
      type(key_set_t), intent(inout) :: keys
      character(len=*), intent(in) :: name

      if (keys%n_required >= MAX_KEYS) then
         error stop "mqc_json_schema: required-list is full; raise MAX_KEYS"
      end if
      keys%n_required = keys%n_required + 1
      keys%required(keys%n_required) = name
   end subroutine require

   pure function is_allowed(keys, name) result(ok)
      type(key_set_t), intent(in) :: keys
      character(len=*), intent(in) :: name
      logical :: ok
      integer :: i

      ok = .false.
      do i = 1, keys%n_allowed
         if (trim(keys%allowed(i)) == trim(name)) then
            ok = .true.
            return
         end if
      end do
   end function is_allowed

   pure function allowed_list(keys) result(text)
      !! The allowed keys, for an error message that says what to write instead
      type(key_set_t), intent(in) :: keys
      character(len=:), allocatable :: text
      integer :: i

      text = ""
      do i = 1, keys%n_allowed
         if (i > 1) text = text//", "
         text = text//trim(keys%allowed(i))
      end do
   end function allowed_list

   pure function quoted(path, name) result(text)
      !! A key named as the user would write it, with its parent path
      character(len=*), intent(in) :: path, name
      character(len=:), allocatable :: text

      if (len_trim(path) == 0) then
         text = "'"//trim(name)//"'"
      else
         text = "'"//trim(path)//"."//trim(name)//"'"
      end if
   end function quoted

   pure function molecule_path(imol) result(text)
      character(len=*), parameter :: PREFIX = "molecules["
      integer, intent(in) :: imol
      character(len=:), allocatable :: text

      ! Written 0-based to match how a user counts entries in the array.
      text = PREFIX//int_to_text(imol - 1)//"]"
   end function molecule_path

   function has_key(core, object, name) result(present_key)
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: object
      character(len=*), intent(in) :: name
      logical :: present_key

      type(json_value), pointer :: child

      call core%get(object, name, child, present_key)
      present_key = present_key .and. associated(child)
   end function has_key

   function child_count(core, object, name) result(n)
      !! Entries in `object.name`, or 0 when it is absent
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: object
      character(len=*), intent(in) :: name
      integer :: n

      type(json_value), pointer :: child
      logical :: found

      n = 0
      call core%get(object, name, child, found)
      if (.not. found .or. .not. associated(child)) return
      call core%info(child, n_children=n)
   end function child_count

   function sum_integer_array(core, array, n) result(total)
      type(json_core), intent(inout) :: core
      type(json_value), pointer, intent(in) :: array
      integer, intent(in) :: n
      integer :: total

      type(json_value), pointer :: element
      integer :: i, value
      logical :: found

      total = 0
      do i = 1, n
         call core%get_child(array, i, element, found)
         if (.not. found) cycle
         call core%get(element, value)
         total = total + value
      end do
   end function sum_integer_array

   pure function nmer_level(name) result(level)
      !! The n-mer level a cutoff key names, or 0 if it names none
      character(len=*), intent(in) :: name
      integer :: level

      integer :: status

      select case (trim(name))
      case ("dimer")
         level = 2
      case ("trimer")
         level = 3
      case ("tetramer")
         level = 4
      case ("pentamer")
         level = 5
      case ("hexamer")
         level = 6
      case ("heptamer")
         level = 7
      case ("octamer")
         level = 8
      case default
         read (name, *, iostat=status) level
         if (status /= 0 .or. level < 2) level = 0
      end select
   end function nmer_level

end module mqc_json_schema
