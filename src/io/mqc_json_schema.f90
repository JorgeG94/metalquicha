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

   integer, parameter :: MAX_KEYS = 20
      !! Widest allow-list below; raising it is the only cost of a new key.
      !! fragmentation_keys is the widest at 18 (its FMO inner-SCF controls
      !! pushed it past 16), so this keeps a little headroom above that.

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
      call check_grandchild_object(core, root, "properties", "bonding_analysis", &
                                   bonding_analysis_keys(), error)
      call check_bonding_analysis(core, root, error)
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
      call allow(keys, "aux_basis")
      call allow(keys, "functional")
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
   end function properties_keys

   function fukui_keys() result(keys)
      !! Settings for the Fukui index analysis
      !!
      !! `population` is required rather than defaulted, because the object
      !! being present is what asks for the analysis and because the scheme
      !! changes the numbers: a choice that large should be written down in the
      !! deck rather than inherited silently.
      type(key_set_t) :: keys
      call allow(keys, "population")
      call require(keys, "population")
   end function fukui_keys

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
      call allow(keys, "unrestricted")
      call allow(keys, "guess")
      call allow(keys, "allow_crap_scf")
      call allow(keys, "density_fitting")
   end function scf_keys

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
      call allow(keys, "radial_points")
      call allow(keys, "angular_points")
   end function dft_keys

   function pcm_keys() result(keys)
      !! The polarizable continuum: the cavity, the solvent, and the charge solve
      !!
      !! Separate from `xtb`, which configures tblite's own CPCM. Two continuum
      !! implementations with different cavities should not share one keyword.
      type(key_set_t) :: keys
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
                        int_text(n_spaces)//" subspaces but gives "// &
                        int_text(n_min)//" minima and "//int_text(n_max)// &
                        " maxima. All three lists describe the same subspaces and "// &
                        "have to be the same length.")
         return
      end if
   end subroutine check_ormas_group

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
      type(key_set_t), intent(inout) :: keys
      character(len=*), intent(in) :: name

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
